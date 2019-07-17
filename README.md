# GeneShot

A gene-level metagenomics oriented pipeline for the analysis of shotgun (WGS) reads from microbial communities.

This repository contains two pipelines, Sciluigi and Nextflow. Both of these workflow management
systems are intended to run software seamlessly on local computer resources (via docker), on
private HPC clusters (using singularity), or on public clouds (e.g. AWS batch).

## Running Nextflow Pipeline

The analysis pipeline can be run with a single command:

```
# Reference database
REF_DMND="s3://fh-ctr-public-reference-data/tool_specific_data/cagalog/2019-06-26-CAGalog/2019-06-26-CAGalog.dmnd"
REF_HDF="s3://fh-ctr-public-reference-data/tool_specific_data/cagalog/2019-06-26-CAGalog/2019-06-26-CAGalog.hdf5"

nextflow \
    -C ~/nextflow.config \
    run \
    jgolob/geneshot \
        --batchfile $BATCHFILE \
        --paired \
        --ref_dmnd $REF_DMND \
        --ref_hdf5 $REF_HDF5 \
        --output_folder $OUTPUT_FOLDER \
        --output_prefix $PROJECT \
        -with-report $PROJECT.html \
        -work-dir $WORK_DIR \
        -resume
```

### Nextflow Configuration File

The argument `~/nextflow.config` refers to the 
[Nextflow configuration file](https://www.nextflow.io/docs/latest/config.html) 
which is typically created when first 
[installing and setting up Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation). 
This is the file that you will use to set up 
[execution on AWS](https://www.nextflow.io/docs/latest/config.html#config-aws), for example.


### Batchfile

The location of all input data is given with the `--batchfile` parameter. Depending
on how the input data is formatted, that file will specify either the single-end, paired-end,
or interleaved FASTQ data. Alternatively, you can also specify a set of SRA accessions, which
will be downloaded and analyzed. 

The default setting is for single-end FASTQ data, in which case the batch file will contain
at least two columns: `name` and `fastq`. Including additional columns is not a problem.

In the example above, the `--paired` flag is set, which indicates that the batchfile
will contain three columns: `name`, `fastq1`, and `fastq2`.

For interleaved data, use the `--interleaved` flag, and use the columns `name` and `fastq`.

For data in NCBI's SRA repository, use the `--from_ncbi_sra` flag, and include the columns
`name` and `run` (with the SRA accession).

The batchfile will be included in the final output file, so it's a convenient place to store
metadata describing each sample.


### Reference Database

The reference database consists of two objects, a database of microbial genes which has been
indexed for alignment by DIAMOND, and an HDF5 store containing metadata describing those genes.

These files can be accessed at the following bucket, or just referenced directly at runtime:

  * `s3://fh-ctr-public-reference-data/tool_specific_data/cagalog/2019-06-26-CAGalog/2019-06-26-CAGalog.dmnd`
  * `s3://fh-ctr-public-reference-data/tool_specific_data/cagalog/2019-06-26-CAGalog/2019-06-26-CAGalog.hdf5`


### Output Folder

The `--output_folder` flag specifies a folder in which all output files will be placed.
The `--output_prefix` flag specifies the prefix that will be appended to all output files.


### HUMAnN2

To run HUMAnN2 on the same datasets, use the `--humann` flag.


### Execution Report

Use the `-with-report` flag to write out a summary of the analysis execution as HTML.


### Working Directory

Use the `-work-dir` flag to specify the working directory for Nextflow execution 
([docs](https://www.nextflow.io/docs/latest/index.html)).


### Call Caching

Use the `-resume` flag to [automatically retrieve cached calls](https://www.nextflow.io/docs/latest/process.html#cache).


## Output Data Format

The output files from geneshot is provided in two formats. Each set of outputs is provided
both as CSV and in a single HDF file. The HDF file is the canonical output format, as it 
contains all of the output data in a single file. Below are the objects that can be found 
within the HDF store.


#### Metadata


Key: `/metadata`

Description: The contents of the batchfile provided by the user.

Example:

| run | name | Cancer |
|  -----  |  -----  |  -----  |
| ERR1018185 | ERR1018185 | 1 |
| ERR1018186 | ERR1018186 | 1 |
| ERR1018187 | ERR1018187 | 1 |
| ERR1018188 | ERR1018188 | 1 |
| ERR1018189 | ERR1018189 | 1 |

Key: `/readcounts`

Description: The number of reads, both total and aligned, for each sample

Example:

| name | n_reads | aligned_reads |
|  -----  |  -----  |  -----  |
| ERR1018185 | 60749206 | 20735371 |
| ERR1018186 | 52854792 | 7869092 |
| ERR1018188 | 55761702 | 19732306 |
| ERR1018187 | 59873382 | 21136212 |
| ERR1018190 | 52102228 | 4358801 |

Key: `/abund/alleles`

Description: The abundance of each microbial gene (`id`) across each sample, including the
coverage, depth of sequencing, length of the gene, number of reads aligned to each gene, 
as well as the NCBI taxid for each gene.

Example:

| coverage | depth | id | length | nreads | std | sample | prop | taxid |
|  -----  |  -----  |  -----  |  -----  |  -----  |  -----  |  -----  |  -----  |  -----  |
| 0.4782608695652174 | 0.4782608695652174 | MGYA00137305.46836 | 69 | 1 | 0.4995271866554809 | ERR1018287 | 1.2715449280078178e-07 | 186802 |
| 0.7402597402597403 | 16.01948051948052 | MGYA00241128.12019 | 154 | 84 | 13.456264306292045 | ERR1018287 | 4.259075015353695e-06 | 0 |
| 1.0 | 5.03125 | MGYA00241455.14143 | 64 | 12 | 1.530612112032307 | ERR1018287 | 1.3376508148900423e-06 | 29347 |
| 0.9485981308411215 | 3.485981308411215 | MGYA00198710.15018 | 214 | 26 | 2.0110186472657645 | ERR1018287 | 9.268125690211613e-07 | 292800 |
| 0.5730337078651685 | 2.067415730337079 | MGYA00140016.37246 | 89 | 6 | 2.1242196336272263 | ERR1018287 | 5.49660687980499e-07 | 39491 |


Key: `/abund/CAGs`

Description: The proportional abundance of each CAG in each sample. To make this calculation, 
the depth of reads aligned to each `allele` are summed according to the `gene` they belng to, 
and then that value is averaged across all `genes` in the same CAG (`group`, below). 

Example:

| sample | group | prop |
|  -----  |  -----  |  -----  |
| ERR1018185 | 0.0 | 4.4889710716337803e-07 |
| ERR1018185 | 1.0 | 2.4594611600221414e-07 |
| ERR1018185 | 2.0 | 3.466341700459156e-07 |
| ERR1018185 | 3.0 | 3.755588272415972e-07 |
| ERR1018185 | 4.0 | 2.0103208922176596e-07 |


Key: `/abund/KEGG_KO`

Description: The proportional abundance of all genes with a given KEGG annotation, 
in each sample.

Example:

| sample | KO | prop | nreads |
|  -----  |  -----  |  -----  |  -----  |
| ERR1018185 | K00001 | 0.00016733925639759165 | 2561 |
| ERR1018185 | K00002 | 8.611910744418914e-06 | 117 |
| ERR1018185 | K00003 | 0.00015178288314809756 | 2499 |
| ERR1018185 | K00004 | 0.0002231820584130467 | 3744 |
| ERR1018185 | K00005 | 5.8290034832123886e-05 | 846 |


Key: `/abund/metaphlan`

Description: The abundance of microbial taxa in each sample, as assayed by MetaPhlAn2.

Example:

| abund | rank | org_name | sample |
|  -----  |  -----  |  -----  |  -----  |
| 1.0 | kingdom | Bacteria | ERR1018204 |
| 0.6298709 | phylum | Verrucomicrobia | ERR1018204 |
| 0.29485190000000006 | phylum | Firmicutes | ERR1018204 |
| 0.0314851 | phylum | Bacteroidetes | ERR1018204 |
| 0.0309732 | phylum | Proteobacteria | ERR1018204 |


Key: `/annot/alleles`

Description: The annotation of each individual `allele`, both KEGG KO and NCBI taxid.

Example:

| allele | KOs | taxid |
|  -----  |  -----  |  -----  |
| MGYA00002187.1293 | nan | 0 |
| MGYA00002318.11 | nan | 0 |
| MGYA00003148.1507 | nan | 0 |
| MGYA00012323.10 | nan | 68766 |
| MGYA00012323.33 | nan | 848 |


Key: `/groups/CAGs`

Description: Groupings of alleles into genes and CAGs.

Example:

| allele | gene | group |
|  -----  |  -----  |  -----  |
| MGYA00137475.60798 | MGYA00137475.60798 | 49088 |
| MGYA00154982.53443 | MGYA00137475.60798 | 49088 |
| MGYA00139858.21391 | MGYA00139858.21391 | 79131 |
| MGYA00242760.39223 | MGYA00242760.39223 | 6033 |
| MGYA00243832.15095 | MGYA00242760.39223 | 6033 |


Key: `/groups/KEGG_KO`

Description: The KEGG KO annotation for each allele.

Example:

| KO | allele |
|  -----  |  -----  |
| K01153 | MGYA00065482.1571 |
| K00558 | MGYA00065482.15104 |
| K03607 | MGYA00065483.5649 |
| K03657 | MGYA00065483.5692 |
| K02016 | MGYA00065483.13065 |


Key: `/groups/NCBI_TAXID`

Description: The NCBI taxid annotation for each allele.

Example:

| allele | taxid |
|  -----  |  -----  |
| MGYA00065482.1541 | 0 |
| MGYA00065482.1545 | 0 |
| MGYA00065482.1563 | 0 |
| MGYA00065482.1568 | 0 |
| MGYA00065482.1570 | 1301 |

