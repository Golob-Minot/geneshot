#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl
import os


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

class RemoveAdapters(sl.ContainerTask):
    container = 'golob/cutadapt:1.18__bcw.0.3.0_al38'
    in_reads = None
    trimmed_R1_path = sl.Parameter()
    trimmed_R2_path = sl.Parameter()
    cutadapt_log_path = sl.Parameter()
    # Default adapter sequences are for nextera, and from:
    # https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
    adapter_F = sl.Parameter(default='CTGTCTCTTATACACATCT')
    adapter_R = sl.Parameter(default='AGATGTGTATAAGAGACAG')

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
    container = 'golob/spades:3.12.0__bcw.0.3.0'

    in_reads = None

    destination_dir = sl.Parameter()
    container_temp_dir = sl.Parameter(default='/scratch/')

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
                'mkdir -p /working &&'
                'metaspades.py '
                '--meta '
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
                'mem': self.containerinfo.mem / 1024,
            }
        )


# Workflow
class Workflow_SGOM(sl.WorkflowTask):
    def workflow(self):
        pass
        # We need to load in a manifest, a CSV file with a line per pair of reads. 
        # A given specimen (a microbial community) can have MULTIPLE paired reads.
        # For each pair of reads, we need to:
        # - QC with cutadapt to remove adapter seqs and low-quality reads
        # - (Optionally): Remove human reads

        # Then *by specimen* we need to:
        # - combine reads into a single pair of
        # forward and revese files for the specimen
        # - Assemble (metaspades)
        # - extract 16S (emirge)
        # - Compositional determination (Metaphlan2 / kraken / etc)


class SHOTGUNOMATIC:
    def __init__(self):
        # Make a parser
        parser = argparse.ArgumentParser(description="""
        Basic Pipeline for processing reads from shotgun-microbiome""")

        if len(sys.argv) < 2:
            parser.print_help()
        else:
            # Common options here
            parser.add_argument(
                '--luigi-manager',
                help="""Run with luigi's work manager daemon (default False)""",
                action='store_true',
            )

            parser.add_argument(
                '--workers',
                help="""How many concurrent workers to use""",
                type=int,
                default=1,
            )
