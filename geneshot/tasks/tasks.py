import logging
import sys
import luigi
import sciluigi as sl
import os
import csv
from collections import defaultdict
import uuid
import gzip
import tarfile

log = logging.getLogger('sciluigi-interface')


# Tasks
class LoadFile(sl.ExternalTask):
    path = sl.Parameter()
    file_format = sl.Parameter(default=None)

    def out_file(self):
        if self.file_format == 'gzip':
            file_format = luigi.format.Gzip
        elif self.file_format == 'Nop':
            file_format = luigi.format.Nop
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

class BCCSpecimenReads(sl.ContainerTask):
    # Use barcodecop to verify reads were properly demultiplexed
    container = 'golob/barcodecop:0.4.1__bcw_0.3.0'

    # For dependencies
    in_reads = None

    specimen = sl.Parameter()
    path = sl.Parameter()

    def out_reads(self):
        reads_dict = {}
        reads_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.bcc.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        reads_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.bcc.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return reads_dict
    """
    def out_stats(self):
        stats_dict = {}
        stats_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.bcc.stats.csv".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        stats_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.bcc.stats.csv".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return stats_dict
    """

    def run(self):
        input_targets={
            'read_1': self.in_reads().get('R1'),
            'read_2': self.in_reads().get('R2'),
            'index_1': self.in_reads().get('I1'),
            'index_2': self.in_reads().get('I2'),
        }
        output_targets={
            'bcc_read_1': self.out_reads()['R1'],
            'bcc_read_2': self.out_reads()['R2'],
            #'bcc_stat_1': self.out_stats()['R1'],
            #'bcc_stat_2': self.out_stats()['R2']
        }
        self.ex(
            command=(
                'barcodecop'
                ' $index_1 $index_2'
                ' --match-filter'
                ' -f $read_1'
                ' -o $bcc_read_1'
            #    ' -C $bcc_stat_1'
                ' && barcodecop'
                ' $index_1 $index_2'
                ' --match-filter'
                ' -f $read_2'
                ' -o $bcc_read_2'
            #    ' -C $bcc_stat_2'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )

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


class BWAIndexHumanGenome(sl.ContainerTask):
    container = 'golob/bwa:0.7.17__bcw.0.3.0C'

    genome_source = sl.Parameter(
        default='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz'
    )
    genome_index_source = sl.Parameter(
        default='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
    )
    exclude_EBV = sl.Parameter(
        default=False
    )
    exclude_MT = sl.Parameter(
        default=False
    )
    index_tgz_path = sl.Parameter()

    def out_bwa_index(self):
        return sl.ContainerTargetInfo(
            self,
            self.index_tgz_path,
            format=luigi.format.Nop
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
        output_targets = {
            'index_tgz': self.out_bwa_index()
        }
        if self.exclude_EBV is False and self.exclude_MT is False:
            # We are just going to keep the complete genome build
            self.ex(
                command=(
                    'wget $genome_index_source -O $index_tgz '
                ),
                output_targets=output_targets,
                extra_params={
                    'genome_index_source': self.genome_index_source
                }
            )
        else:  # We ARE excluding portions of the genome
            exclusion_str = ""
            if self.exclude_EBV and self.exclude_MT:
                exclusion_str = "if sr.id == 'chrEBV' or sr.id == 'chrM': continue"
            elif self.exclude_EBV:
                exclusion_str = "if sr.id == 'chrEBV': continue"
            elif self.exclude_MT:
                exclusion_str = "if sr.id == 'chrM': continue"

            self.ex(
                command=(
                    'mkdir -p $working_dir && '
                    'mkdir -p $working_dir/bwa_index &&'
                    'cd $working_dir && '
                    'wget $genome_source -O $working_dir/genome.fna.gz && '
                    'echo -e "'
                    'from Bio import SeqIO\n'
                    'import gzip\n'
                    "ofh = gzip.open('$working_dir/genome.filtered.fna.gz', 'wt')\n"
                    "for sr in SeqIO.parse(gzip.open('$working_dir/genome.fna.gz', 'rt'), 'fasta'):\n"
                    "\t$exclusion_str\n"
                    "\tSeqIO.write(sr, ofh, 'fasta')\n"
                    "\n"
                    "ofh.close()\n"
                    '" | python && '
                    'cd $working_dir/bwa_index && '
                    'bwa index -a bwtsw -p human_genome_bwa $working_dir/genome.filtered.fna.gz && '
                    'tar czvf $index_tgz . '
                ),
                output_targets=output_targets,
                extra_params={
                    'working_dir': container_working_dir,
                    'genome_source': self.genome_source,
                    'exclusion_str': exclusion_str,
                }
            )


class RemoveHuman(sl.ContainerTask):
    # Human genome build can be found at:
    # ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/

    container = 'golob/bwa:0.7.17__bcw.0.3.0C'
    # Dependencies / inputs
    in_reads = None
    in_bwa_index = None

    # Path parameters
    R1_path = sl.Parameter()
    R2_path = sl.Parameter()
    log_path = sl.Parameter()

    # Alignment parameters
    min_align_score = sl.Parameter(default=30)

    def out_reads(self):
        return {
            'R1': sl.ContainerTargetInfo(
                self,
                self.R1_path,
                format=luigi.format.Nop
            ),
            'R2': sl.ContainerTargetInfo(
                self,
                self.R2_path,
                format=luigi.format.Nop
            ),
        }

    def out_log(self):
        return sl.ContainerTargetInfo(
            self,
            self.log_path,
            format=luigi.format.Nop
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2'],
            'bwa_index': self.in_bwa_index()
        }
        output_targets = {
            'unaligned_1': self.out_reads()['R1'],
            'unaligned_2': self.out_reads()['R2'],
            'log': self.out_log()
        }
        with self.in_bwa_index().open('rb') as bwa_i_h:
            bwa_i_tar = tarfile.open(
                fileobj=bwa_i_h,
                mode='r'
            )
            first_in_index = bwa_i_tar.next()
            first_in_index_sp = first_in_index.name.split('.')
            bwa_index_prefix = ".".join(first_in_index_sp[0:len(first_in_index_sp)-1])

        self.ex(
            command=(
                'mkdir -p $working_dir && '
                'cd $working_dir && '
                'echo BWA index file prefix is $bwa_index_prefix | tee -a $log && '
                'tar xzvf $bwa_index | tee - a $log && '
                'echo Files in working directory: | tee - a $log && '
                'ls -lh | tee - a $log && '
                'echo Running BWA | tee - a $log && '
                'bwa mem '
                '-t $vCPU '
                '-o alignment.sam '
                '$bwa_index_prefix '
                '$read_1 '
                '$read_2 | tee - a $log && '
                'echo Extracting Unaligned Pairs && '
                'samtools fastq alignment.sam -f 12 --threads $vCPU -1 $unaligned_1 -2 $unaligned_2 > $log && '
                'echo Deleting temporary folder | tee -a $log && '
                'rm -r $working_dir && '
                'echo Done | tee -a $log'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': self.containerinfo.vcpu,
                'working_dir': container_working_dir,
                'min_aln_score': self.min_align_score,
                'bwa_index_prefix': bwa_index_prefix
            }
        )


class AlignReads(sl.ContainerTask):
    container = 'golob/bwa:0.7.17__bcw.0.3.0C'
    in_reads = None
    # Human genome build can be found at:
    # ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
    in_bwa_index = None
    bam_path = sl.Parameter()
    alignment_log_path = sl.Parameter()

    def out_bam(self):
        return sl.ContainerTargetInfo(
            self,
            self.bam_path,
            format=luigi.format.Nop
        )

    def out_log(self):
        return sl.ContainerTargetInfo(
            self,
            self.alignment_log_path,
            format=luigi.format.Nop
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
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
                'bwa_index_prefix=$$(tar -ztvf $bwa_index | head -1 | sed \'s/.* //\' | sed \'s/.amb//\') && '
                'echo BWA index file prefix is $${bwa_index_prefix} | tee - a $log && '
                'tar xzvf $bwa_index | tee - a $log && '
                'echo Files in working directory: | tee - a $log && '
                'ls -lh | tee - a $log && '
                'echo Running BWA | tee - a $log && '
                'bwa mem '
                '-t $vCPU '
                '-o alignment.sam '
                '$${bwa_index_prefix} '
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
                'working_dir': container_working_dir,
            }
        )


class ExtractUnalignedPairs(sl.ContainerTask):
    container = 'golob/bwa:0.7.17__bcw.0.3.0C'
    in_bam = None
    unaligned_R1_path = sl.Parameter()
    unaligned_R2_path = sl.Parameter()
    unaligned_log_path = sl.Parameter()

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
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
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
                'samtools fastq $bam -f 12 --threads $vCPU -1 $unaligned_1 -2 $unaligned_2 > $log'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vCPU': str(int(self.containerinfo.vcpu) - 1),
                'working_dir': container_working_dir,
            }
        )


class CombineReads(sl.ContainerTask):
    container = 'golob/fastatools:0.6.2__bcw.0.3.0'

    in_reads_list = None
    combined_R1_path = sl.Parameter()
    combined_R2_path = sl.Parameter()

    def out_reads(self):
        return {
            'R1': sl.ContainerTargetInfo(
                self,
                self.combined_R1_path,
                format=luigi.format.Nop
            ),
            'R2': sl.ContainerTargetInfo(
                self,
                self.combined_R2_path,
                format=luigi.format.Nop
            ),
        }

    def run(self):
        assert len(self.in_reads_list) >= 2, "Less than two read pairs to combine. Nothing to do"
        # Implicit else
        input_targets = {}
        R1_command_string = ""
        R2_command_string = ""
        for index, read_pair in enumerate(self.in_reads_list):
            input_targets['p{}_R1'.format(index)] = read_pair()["R1"]
            input_targets['p{}_R2'.format(index)] = read_pair()["R2"]
            R1_command_string += '$p{}_R1 '.format(index)
            R2_command_string += '$p{}_R2 '.format(index)
        command = (
            'combine_fastq_pairs.py '
            '-1 {} '
            '-2 {} '
            '--normalize-ids '
            '-o1 $out_R1 '
            '-o2 $out_R2'
        ).format(
            R1_command_string,
            R2_command_string
        )
        self.ex(
            command=command,
            input_targets=input_targets,
            output_targets={
                'out_R1': self.out_reads()['R1'],
                'out_R2': self.out_reads()['R2']
            },
            extra_params={
                'R1_string': R1_command_string,
                'R2_string': R2_command_string,
            }
        )


class MetaSPAdesAssembly(sl.ContainerTask):
    # runMetaSpades

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/spades:3.12.0__bcw.0.3.0A'

    in_reads = None

    destination_dir = sl.Parameter()

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
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
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
                'working_dir': container_working_dir,
                'phred_offset': self.phred_offset
            }
        )


class EggnogMapperMap(sl.ContainerTask):
    container = 'golob/eggnog-mapper:1.0.3__bcw.0.3.1'

    # Where to put the gzipped annotations
    annotation_path_gz = sl.Parameter()
    in_db_tgz = None
    in_faa = None

    def out_annotations(self):
        return sl.ContainerTargetInfo(
            self,
            self.annotation_path_gz,
            format=luigi.format.Nop
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
        self.ex(
            command=(
                'mkdir -p $working_dir/emdb && mkdir -p $working_dir/egm/ '
                '&& tar -C $working_dir/emdb/ -xzf $db_tgz ./eggnog.db ./eggnog_proteins.dmnd '
                '&& emapper.py '
                '-i $faa '
                '-m diamond '
                '--dmnd_db $working_dir/emdb/eggnog_proteins.dmnd '
                '--data_dir $working_dir/emdb/ '
                '--cpu $vcpu '
                '--temp_dir $working_dir '
                '-o $working_dir/egm '
                '&& gzip -c $working_dir/egm.emapper.annotations > $annot_gz '
                '; rm -r $working_dir'
            ),
            input_targets={
                'faa': self.in_faa(),
                'db_tgz': self.in_db_tgz(),
            },
            output_targets={
                'annot_gz': self.out_annotations(),
            },
            extra_params={
                'vcpu': self.containerinfo.vcpu,
                'working_dir': container_working_dir,
            }
        )


class EggnogMapperDownloadDB(sl.ContainerTask):
    container = 'golob/eggnog-mapper:1.0.3__bcw.0.3.0'
    db = sl.Parameter(default='none')

    destination_tgz = sl.Parameter()

    def out_eggnog_db_tgz(self):
        return sl.ContainerTargetInfo(
            self,
            self.destination_tgz,
            format=luigi.format.Nop
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
        self.ex(
            command=(
                'mkdir -p $working_dir && '
                'download_eggnog_data.py $db -y --data_dir $working_dir && '
                'tar czvf $eggnog_db_tgz --directory $working_dir . && '
                'rm -r $working_dir'
            ),
            output_targets={
                'eggnog_db_tgz': self.out_eggnog_db_tgz()
            },
            extra_params={
                'working_dir': container_working_dir,
                'db': self.db
            }
        )


class ProkkaAnnotate(sl.ContainerTask):
    container = 'golob/metaspades:v3.11.1--8B__bcw__0.3.0'
    in_contigs = None
    destination_dir = sl.Parameter()
    prefix = sl.Parameter(default='prokka')
    center = sl.Parameter(default='geneshot')

    def out_tbl(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.tbl.gz".format(self.prefix)
            )
        )

    def out_err(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.err.gz".format(self.prefix)
            )
        )

    def out_faa(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.faa.gz".format(self.prefix)
            )
        )

    def out_ffn(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.ffn.gz".format(self.prefix)
            )
        )

    def out_fna(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.fna.gz".format(self.prefix)
            )
        )

    def out_fsa(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.fsa.gz".format(self.prefix)
            )
        )

    def out_gbk(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.gbk.gz".format(self.prefix)
            )
        )

    def out_gff(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.gff.gz".format(self.prefix)
            )
        )

    def out_log(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.log.gz".format(self.prefix)
            )
        )

    def out_sqn(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.sqn.gz".format(self.prefix)
            )
        )

    def out_tsv(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.tsv.gz".format(self.prefix)
            )
        )

    def out_txt(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.destination_dir,
                "{}.txt.gz".format(self.prefix)
            )
        )

    def run(self):
        container_working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4().hex)
        )
        self.ex(
            command=(
                'mkdir -p $container_working_dir '
                '&& ls -l $contigs '
                '&& prokka '
                '--outdir $container_working_dir '
                '--centre $center '
                '--compliant '
                '--prefix prokka '
                '--cpus $vcpu '
                '--force '
                '$contigs '
                '&& gzip -c $container_working_dir/prokka.tbl > $tbl '
                '&& gzip -c $container_working_dir/prokka.err > $err '
                '&& gzip -c $container_working_dir/prokka.faa > $faa '
                '&& gzip -c $container_working_dir/prokka.ffn > $ffn '
                '&& gzip -c $container_working_dir/prokka.fna > $fna '
                '&& gzip -c $container_working_dir/prokka.fsa > $fsa '
                '&& gzip -c $container_working_dir/prokka.gbk > $gbk '
                '&& gzip -c $container_working_dir/prokka.gff > $gff '
                '&& gzip -c $container_working_dir/prokka.log > $log '
                '&& gzip -c $container_working_dir/prokka.sqn > $sqn '
                '&& gzip -c $container_working_dir/prokka.tsv > $tsv '
                '&& gzip -c $container_working_dir/prokka.txt > $txt '
                '&& rm -r $container_working_dir'
            ),
            input_targets={
                'contigs': self.in_contigs(),
            },
            output_targets={
                'tbl': self.out_tbl(),
                'err': self.out_err(),
                'faa': self.out_faa(),
                'ffn': self.out_ffn(),
                'fna': self.out_fna(),
                'fsa': self.out_fsa(),
                'gbk': self.out_gbk(),
                'gff': self.out_gff(),
                'log': self.out_log(),
                'sqn': self.out_sqn(),
                'tsv': self.out_tsv(),
                'txt': self.out_txt(),
            },
            extra_params={
                'vcpu': self.containerinfo.vcpu,
                'container_working_dir': container_working_dir,
                'center': self.center,
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
