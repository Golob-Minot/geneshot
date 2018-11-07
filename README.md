# geneshot

A gene-level metagenomics oriented pipeline for the analysis of shotgun reads from microbial communities.

This pipeline uses sciluigi (customized fork running tasks in containers) to run software seamlessly 
on local computer resources (via docker), on private HPC clusters (using singularity), or on public clouds 
(e.g. AWS batch).

### Config File

All of the parameters needed to run jobs on Slurm, PBS, or AWS Batch are saved in a config file
that contains account information specific to your infrastructure. An example of the sciluigi
config can be found in `geneshot/sciluigi.ini`. Copy that to `~/.sciluigi/` and you can use it
across multiple projects.
