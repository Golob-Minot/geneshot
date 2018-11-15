import logging
import sys
import luigi
import sciluigi as sl
import os
import csv
from collections import defaultdict
import uuid

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


class AlignReads(sl.ContainerTask):
    container = 'quay.io/fhcrc-microbiome/bwa:v0.7.17--4'
    in_reads = None
    in_bwa_index = None
    bam_path = sl.Parameter()
    alignment_log_path = sl.Parameter()
    container_working_dir = sl.Parameter(default=os.path.join(
        '/tmp',
        str(uuid.uuid4())
    ))

    def out_bam(self):
        return sl.ContainerTargetInfo(
                self,
                self.bam_path,
                format=luigi.format.Nop
            )

    def out_log(self):
        return sl.ContainerTargetInfo(
                self,
                self.bam_path,
                format=luigi.format.Nop
            )

    def run(self):
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2'],
            'bwa_index': self.in_bwa_index()
        }
        output_targets = {
            'bam': self.out_bam(),
            'log': self.out_log()
        }
        self.ex(
            command=(
                'mkdir -p $working_dir && '
                'cd $working_dir && '
                'bwa_index_prefix=$(tar -ztvf $bwa_index | head -1 | sed \'s/.* //\' | sed \'s/.amb//\') && '
                'echo BWA index file prefix is ${bwa_index_prefix} | tee - a $log && '
                'tar xzvf $bwa_index | tee - a $log && '
                'echo Files in working directory: | tee - a $log && '
                'ls -lh | tee - a $log && '
                'echo Running BWA | tee - a $log && '
                'bwa mem '
                '-t $vCPU '
                '-o alignment.sam '
                '${bwa_index_prefix} '
                '$read_1 '
                '$read_2 | tee - a $log && '
                'echo Converting to BAM && '
                'samtools view -bh -o $bam alignment.sam && '
                'echo Deleting temporary folder | tee -a $log && '
                'rm -r $working_dir && '
                'echo Done | tee -a $log'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'working_dir': self.container_working_dir,
            }
        )


class ExtractUnalignedPairs(sl.ContainerTask):
    container = 'quay.io/fhcrc-microbiome/bwa:v0.7.17--4'
    in_bam = None
    unaligned_R1_path = sl.Parameter()
    unaligned_R2_path = sl.Parameter()
    unaligned_log_path = sl.Parameter()
    container_working_dir = sl.Parameter(default=os.path.join(
        '/tmp',
        str(uuid.uuid4())
    ))

    def out_reads(self):
        return {
            'R1': sl.ContainerTargetInfo(
                self,
                self.unaligned_R1_path,
                format=luigi.format.Nop
            ),
            'R2': sl.ContainerTargetInfo(
                self,
                self.unaligned_R2_path,
                format=luigi.format.Nop
            ),
            'log': sl.ContainerTargetInfo(
                self,
                self.unaligned_log_path,
                format=luigi.format.Nop
            )
        }

    def run(self):
        input_targets = {
            'bam': self.in_bam()
        }
        output_targets = {
            'unaligned_1': self.out_reads()['R1'],
            'unaligned_2': self.out_reads()['R2'],
            'log': self.out_reads()['log']
        }
        self.ex(
            command=(
                'mkdir -p $working_dir && '
                'cd $working_dir && '
                'echo Extracting unaligned read pairs | tee -a $log && '
                'samtools view -f 12 $bam > unaligned.sam && '
                'echo Extracting R1 | tee -a $log && '
                'samtools view -f 64 unaligned.sam | samtools fastq - | gzip > $unaligned_1 && '
                'echo There are $(cat $unaligned_1 | wc -l) lines in the unmapped R1 file | tee -a $log && '
                'echo Extracting R2 | tee -a $log && '
                'samtools view -F 64 unaligned.sam | samtools fastq - | gzip > $unaligned_2 && '
                'echo There are $(cat $unaligned_2 | wc -l) lines in the unmapped R2 file | tee -a $log && '
                'echo Deleting temporary folder && '
                'rm -r $working_dir && '
                'echo Done | tee -a $log'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'working_dir': self.container_working_dir,
            }
        )


class MetaSPAdesAssembly(sl.ContainerTask):
    # runMetaSpades

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/spades:3.12.0__bcw.0.3.0A'

    in_reads = None

    destination_dir = sl.Parameter()
    container_working_dir = sl.Parameter(default=os.path.join(
        '/tmp',
        str(uuid.uuid4())
    ))
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
                'mkdir -p $working_dir && '
                'metaspades.py '
                '--meta '
                '--phred-offset $phred_offset '
                '-1 $read_1 '
                '-2 $read_2 '
                '-o $working_dir '
                '-t $vCPU '
                '--memory $mem '
                '&& cp $working_dir/contigs.fasta $contigs '
                '&& cp $working_dir/scaffolds.fasta $scaffolds '
                '&& rm -r $working_dir'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'mem': int(self.containerinfo.mem / 1024),
                'working_dir': self.container_working_dir,
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
