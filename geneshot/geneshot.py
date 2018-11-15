#!/usr/bin/env python3
import argparse
import logging
import sys
import luigi
import sciluigi as sl
import os
import csv
from collections import defaultdict
from tasks.tasks import LoadFile
from tasks.tasks import LoadManifest
from tasks.tasks import LoadPairedReads
from tasks.tasks import RemoveAdapters
from tasks.tasks import ExtractUnalignedPairs
from tasks.tasks import AlignReads
from tasks.tasks import MetaSPAdesAssembly
from tasks.tasks import ProkkaAnnotate
from tasks.tasks import EggnogMapperDownloadDB

log = logging.getLogger('sciluigi-interface')

# Workflow
class Workflow_SGOM(sl.WorkflowTask):
    slconfig = sl.Parameter()
    manifest_path = sl.Parameter()
    working_dir = sl.Parameter()
    human_bwa_index = sl.Parameter()

    def workflow(self):
        # Initialize our containerinfo classes from our SL-config file
        light_containerinfo = sl.ContainerInfo(
            # Format is {'source_path': {'bind': '/container/path', 'mode': mode}} 
            mounts={"docker_scratch": {"bind": "/tmp/", "mode": "rw"}}
        )
        light_containerinfo.from_config(
            self.slconfig,
            'light'
        )
<<<<<<< HEAD
        mid_cpu = sl.ContainerInfo()
        mid_cpu.from_config(
            self.slconfig,
            'midcpu'
        )
        heavy_containerinfo = sl.ContainerInfo()
=======
        heavy_containerinfo = sl.ContainerInfo(
            # Format is {'source_path': {'bind': '/container/path', 'mode': mode}}
            mounts={"docker_scratch": {"bind": "/tmp/", "mode": "rw"}}
        )
>>>>>>> align-to-human
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

        # Load the file with the human genome
        human_bwa_index = self.new_task(
            'load_human_bwa_index',
            LoadFile,
            path=self.human_bwa_index,
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

                # Specify the folder for the alignment against the human genome
                sp_read_path_base = os.path.join(
                    self.working_dir,
                    'qc',
                    'noadapt_align_human',
                    '{}.{}'.format(
                        specimen,
                        sp_read_idx
                    )
                )

                # Aligning reads against the human genome, so that the unaligned ones can be removed
                specimen_reads_tasks[specimen]['noadapt_align_human'][sp_read_idx] = self.new_task(
                    'noadapt_align_human.{}.{}'.format(specimen, sp_read_idx),
                    AlignReads,
                    containertargetinfo=heavy_containerinfo,
                    bam_path="{}.bam".format(sp_read_path_base),
                    alignment_log_path="{}.log".format(sp_read_path_base)
                )
                specimen_reads_tasks[specimen]['noadapt_align_human'][sp_read_idx].in_bwa_index = human_bwa_index.out_file
                specimen_reads_tasks[specimen]['noadapt_align_human'][sp_read_idx].in_reads = specimen_reads_tasks[specimen]['noadapt'][sp_read_idx].out_reads

                # Specify the folder for the reads that have had adapters removed, and have also have had human removed                
                sp_read_path_base = os.path.join(
                    self.working_dir,
                    'qc',
                    'noadapt_nohuman',
                    '{}.{}'.format(
                        specimen,
                        sp_read_idx
                    )
                )

                # Aligning reads against the human genome, so that the unaligned ones can be removed
                specimen_reads_tasks[specimen]['noadapt_nohuman'][sp_read_idx] = self.new_task(
                    'noadapt_nohuman.{}.{}'.format(specimen, sp_read_idx),
                    ExtractUnalignedPairs,
                    containertargetinfo=light_containerinfo,
                    unaligned_R1_path="{}.R1.fastq.gz".format(sp_read_path_base),
                    unaligned_R2_path="{}.R2.fastq.gz".format(sp_read_path_base),
                    unaligned_log_path="{}.log".format(sp_read_path_base)
                )
                specimen_reads_tasks[specimen]['noadapt_nohuman'][sp_read_idx].in_bam = specimen_reads_tasks[specimen]['noadapt_align_human'][sp_read_idx].out_bam


            # Combine a given specimen's trimmed and human-depleted reads into one pair of reads TODO
            specimen_combined_reads = specimen_reads_tasks[specimen]['noadapt_nohuman'][0]

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


class GENESHOT:
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

        # Options specific to GeneShot
        parser.add_argument(
            '--manifest',
            '-M',
            help="""Location of manifest file in CSV format""",
            required=True
        )
        # Human genome indexed with BWA and stored in a TGZ
        parser.add_argument(
            '--human-bwa-index',
            help="""Location of human genome BWA index (TGZ)""",
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
            '--human-bwa-index={}'.format(
                args.human_bwa_index
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
    GENESHOT()


if __name__ == "__main__":
    main()
