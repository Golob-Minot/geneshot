#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl
import os
import csv
from collections import defaultdict

log = logging.getLogger('sciluigi-interface')


# Tasks
class LoadFile(sl.ExternalTask):
    path = sl.Parameter()
    file_format = sl.Parameter(default=None)

    def out_file(self):
        if self.file_format == 'gzip':
            file_format = luigi.format.Gzip
        else:
            file_format = None

        return sl.ContainerTargetInfo(self, self.path, format=file_format)


class LoadPairedReads(sl.ExternalTask):
    path_R1 = sl.Parameter()
    path_R2 = sl.Parameter()
    file_format = sl.Parameter(default='gzip')

    def out_reads(self):
        if self.file_format == 'gzip':
            file_format = luigi.format.Gzip
        else:
            file_format = None

        return {
            'R1': sl.ContainerTargetInfo(self, self.path_R1, format=file_format),
            'R2': sl.ContainerTargetInfo(self, self.path_R2, format=file_format)
        }


class LoadManifest(LoadFile):
    def manifest(self):
        # Load manifest as csv and return as a list of dicts
        with self.out_file().open() as manifest_h:
            return [r for r in csv.DictReader(manifest_h)]

    def group_by_specimen(self):
        manifest = self.manifest()
        specimens = {r.get('specimen') for r in manifest}
        log.info("{} specimens".format(
            len(specimens))
        )

        for specimen in specimens:
            yield((specimen, [
                r for r in manifest
                if r.get('specimen') == specimen
            ]))

    def get_columns(self):
        with self.out_file().open() as manifest_h:
            return set(csv.DictReader(manifest_h).fieldnames)

    def is_paired(self):
        return 'R2_path' in self.get_columns()

    def is_valid(self):
        columns = self.get_columns()
        if 'R1_path' in columns and 'specmen' in columns:
            return True
        else:
            return False

    def get_specimens(self):
        with self.out_file().open() as manifest_h:
            return {r.get('specimen') for r in csv.DictReader(manifest_h)}


class RemoveAdapters(sl.ContainerTask):
    container = 'golob/cutadapt:1.18__bcw.0.3.0_al38'
    in_reads = None
    trimmed_R1_path = sl.Parameter()
    trimmed_R2_path = sl.Parameter()
    cutadapt_log_path = sl.Parameter()
    # Default adapter sequences are for nextera, and from:
    # https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
    # Defaulted for nextera
    adapter_F = sl.Parameter(default='CTGTCTCTTATACACATCT')
    adapter_R = sl.Parameter(default='CTGTCTCTTATACACATCT')

    def out_reads(self):
        return {
            'R1': sl.ContainerTargetInfo(
                self,
                self.trimmed_R1_path,
                format=luigi.format.Nop
            ),
            'R2': sl.ContainerTargetInfo(
                self,
                self.trimmed_R2_path,
                format=luigi.format.Nop
            ),
            'log': sl.ContainerTargetInfo(
                self,
                self.cutadapt_log_path,
                format=luigi.format.Nop
            )
        }

    def run(self):
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2']
        }
        output_targets = {
            'trimmed_1': self.out_reads()['R1'],
            'trimmed_2': self.out_reads()['R2'],
            'log': self.out_reads()['log']
        }
        self.ex(
            command=(
                'cutadapt '
                '-j $vCPU '
                '-a $adapter_F '
                '-A $adapter_R '
                '-o $trimmed_1 '
                '-p $trimmed_2 '
                '$read_1 '
                '$read_2 '
                '> $log'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'adapter_F': self.adapter_F,
                'adapter_R': self.adapter_R
            }
        )


class MetaSPAdesAssembly(sl.ContainerTask):
    # runMetaSpades

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/spades:3.12.0__bcw.0.3.0A'

    in_reads = None

    destination_dir = sl.Parameter()
    container_temp_dir = sl.Parameter(default='/scratch/')
    phred_offset = sl.Parameter(default='33')

    def out_contigs(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                'contigs.fasta'
            ),
            format=luigi.format.Nop
        )

    def out_scaffolds(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                'scaffolds.fasta'
            ),
            format=luigi.format.Nop
        )

    def run(self):
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2']
        }
        output_targets = {
            'contigs': self.out_contigs(),
            'scaffolds': self.out_scaffolds(),
        }

        self.ex(
            command=(
                'mkdir -p /working && '
                'metaspades.py '
                '--meta '
                '--phred-offset $phred_offset '
                '-1 $read_1 '
                '-2 $read_2 '
                '-o /working/ '
                '-t $vCPU '
                '--memory $mem '
                '--tmp-dir $tempdir '
                '&& cp /working/contigs.fasta $contigs '
                '&& cp /working/scaffolds.fasta $scaffolds '
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'mem': int(self.containerinfo.mem / 1024),
                'tempdir': self.container_temp_dir,
                'phred_offset': self.phred_offset
            }
        )


# Workflow
class Workflow_SGOM(sl.WorkflowTask):
    slconfig = sl.Parameter()
    manifest_path = sl.Parameter()
    working_dir = sl.Parameter()

    def workflow(self):
        # Initialize our containerinfo classes from our SL-config file
        light_containerinfo = sl.ContainerInfo()
        light_containerinfo.from_config(
            self.slconfig,
            'light'
        )
        heavy_containerinfo = sl.ContainerInfo()
        heavy_containerinfo.from_config(
            self.slconfig,
            'heavy'
        )

        # We need to load in a manifest, a CSV file with a line per pair of reads.
        manifest = self.new_task(
            'load_manifest',
            LoadManifest,
            path=self.manifest_path,
        )

        specimen_reads_tasks = defaultdict(lambda: defaultdict(dict))

        for specimen, specimen_reads in manifest.group_by_specimen():
            # A given specimen (a microbial community) can have MULTIPLE paired reads.
            for sp_read_idx, read_pair in enumerate(specimen_reads):
                # For each pair of reads, we need to:
                # - QC with cutadapt to remove adapter seqs and low-quality reads
                # - (Optionally): Remove human reads
                specimen_reads_tasks[specimen]['raw_reads'][sp_read_idx] = self.new_task(
                    'load.{}.{}'.format(specimen, sp_read_idx),
                    LoadPairedReads,
                    path_R1=read_pair['R1_path'],
                    path_R2=read_pair['R2_path']
                )

                sp_read_path_base = os.path.join(
                    self.working_dir,
                    'qc',
                    'noadapt',
                    '{}.{}'.format(
                        specimen,
                        sp_read_idx
                    )
                )

                specimen_reads_tasks[specimen]['noadapt'][sp_read_idx] = self.new_task(
                    'noadapt.{}.{}'.format(specimen, sp_read_idx),
                    RemoveAdapters,
                    containertargetinfo=light_containerinfo,
                    trimmed_R1_path="{}.R1.fastq.gz".format(sp_read_path_base),
                    trimmed_R2_path="{}.R2.fastq.gz".format(sp_read_path_base),
                    cutadapt_log_path="{}.cutadapt.log".format(sp_read_path_base),
                )
                specimen_reads_tasks[specimen]['noadapt'][sp_read_idx].in_reads = specimen_reads_tasks[specimen]['raw_reads'][sp_read_idx].out_reads
                # remove human here

            # Combine a given specimen's trimmed and human-depleted reads into one pair of reads
            specimen_combined_reads = specimen_reads_tasks[specimen]['noadapt'][0]

            spades_container_info = heavy_containerinfo
            specimen_reads_tasks[specimen]['assembly'] = self.new_task(
                'assemble.{}'.format(specimen),
                MetaSPAdesAssembly,
                containerinfo=heavy_containerinfo,
                destination_dir=os.path.join(
                    self.working_dir,
                    'assembly',
                    specimen
                )
            )

            # - Assemble (metaspades)
            specimen_reads_tasks[specimen]['assembly'].in_reads = specimen_combined_reads.out_reads


            # - extract 16S (emirge)
            # - Compositional determination (Metaphlan2 / kraken / etc)

        return specimen_reads_tasks


class SHOTGUNOMATIC:
    def __init__(self):
        # Make a parser
        parser = argparse.ArgumentParser(description="""
        Basic Pipeline for processing reads from shotgun-microbiome""")

        # Common options here
        parser.add_argument(
            '--luigi-manager',
            help="""Run with luigi's work manager daemon (default False)""",
            action='store_true',
        )

        parser.add_argument(
            '--workers',
            help="""How many concurrent workers to use (default=1)""",
            type=int,
            default=1,
        )

        parser.add_argument(
            '--slconfig',
            help="""Location of sciluigi config file""",
            default=os.path.expanduser('~/.sciluigi/containerinfo.ini')
        )

        # Options specific to shotgunomatic
        parser.add_argument(
            '--manifest',
            '-M',
            help="""Location of manifest file in CSV format""",
            required=True
        )
        parser.add_argument(
            '--working-dir',
            '-w',
            help="""Directory into which we should place our working files""",
            required=True,
        )

        # Unpack CLI options
        args = parser.parse_args()
        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        cmdline_args = [
            '--workers={}'.format(
                args.workers
            ),
            '--slconfig={}'.format(
                args.slconfig
            ),
            '--manifest-path={}'.format(
                args.manifest
            ),
            '--working-dir={}'.format(
                args.working_dir,
            )
        ]

        # Run here
        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=Workflow_SGOM,
            cmdline_args=cmdline_args
        )


def main():
    """Entrypoint for main script."""
    SHOTGUNOMATIC()


if __name__ == "__main__":
    main()
