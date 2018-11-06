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


class Emirge16S(sl.ContainerTask):
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