#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl
import os
import csv
from collections import defaultdict
from tasks.tasks import LoadManifest
from tasks.tasks import LoadPairedReads
from tasks.tasks import RemoveAdapters
from tasks.tasks import MetaSPAdesAssembly
from tasks.tasks import ProkkaAnnotate
from tasks.tasks import EggnogMapperDownloadDB

log = logging.getLogger('sciluigi-interface')

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
        mid_cpu = sl.ContainerInfo()
        mid_cpu.from_config(
            self.slconfig,
            'midcpu'
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

        # download / get eggnog db
        eggnog_dbs = self.new_task(
            'get_eggnog_dbs',
            EggnogMapperDownloadDB,
            destination_tgz=os.path.join(
                self.working_dir,
                'emdb.tgz'
            ),
            containerinfo=light_containerinfo
        )

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

            # - Assemble (metaspades)
            spades_container_info = heavy_containerinfo
            specimen_reads_tasks[specimen]['assembly'] = self.new_task(
                'assemble.{}'.format(specimen),
                MetaSPAdesAssembly,
                containerinfo=spades_container_info,
                destination_dir=os.path.join(
                    self.working_dir,
                    'assembly',
                    specimen
                )
            )
            specimen_reads_tasks[specimen]['assembly'].in_reads = specimen_combined_reads.out_reads

            #  - Annotate Assembly (Prokka)
            specimen_reads_tasks[specimen]['annotate'] = self.new_task(
                'annotate.{}'.format(specimen),
                ProkkaAnnotate,
                containerinfo=mid_cpu,
                destination_dir=os.path.join(
                    self.working_dir,
                    'annotate',
                    'prokka',
                    specimen
                ),
                prefix="".join(s for s in specimen if s.isalnum()),
            )
            specimen_reads_tasks[specimen]['annotate'].in_contigs = specimen_reads_tasks[specimen]['assembly'].out_scaffolds

            #  - eggnog map annotated peptides


            # - extract 16S (emirge)
            # - Compositional determination (Metaphlan2 / kraken / etc)

        return specimen_reads_tasks, eggnog_dbs


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
