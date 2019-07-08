# bismark_snakemake
Snakemake workflow to process and analyse Bisulfite-Sequencing data 

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# ToDo

## Visualization
- visualize with [ViewBS](https://academic.oup.com/bioinformatics/article/34/4/708/4566176) [code](https://github.com/xie186/ViewBS)
- count C coverage aus dem global file, oder den selbsterstellten coverage bedgraph files
- coverage per chromosome (how is the Chloroplast?)
- files mit interessantern Regionen besorgen (gff/gtf, DMRs,  Transposons usw.)

## workflow


### rules
- bismark report rule
- bismark summary report rule
- fastqc rule
- viewBS region rules
- coverage calc rule
- coverage histogram rule
- samtools dict
- dedup rule with https://github.com/FelixKrueger/Bismark/issues/238#issuecomment-470478069
- ViewBS rules output files, not only directories
- multiqc rule

## Snakemake
- define everything  in config file
- conda: for rules (ViewBS for ViewBS, ansonsten ??? gab es schon eines für MS?) dafür müssen die yml files erstellt werden.
- sinnvolle finale all rule
- relative pfade für alle shell commandos
- https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
- local rules localrules: all, clean, make_archive
- scripts have to thow errors/ send correct exit status. for example when output empty or input not correct
- ordnen nach https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling
- apply for https://github.com/snakemake-workflows/

## environment
- rename environment
- create environment file
- if ViewBS does'nt work, remove from environment

## misc
- Documentation
- how to create unit files
