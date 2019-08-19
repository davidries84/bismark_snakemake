
import pandas as pd
configfile: "config.yaml"


SAMPLES = pd.read_table(config["samples"]).set_index("sample", drop=False)
UNITS = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False) # hier werden die Spalten sample und unit im units.tsv als index gesetzt. Muss evtl. erweitert werden
UNITS.index = UNITS.index.set_levels([i.astype(str) for i in UNITS.index.levels]) # indices als str setzten



##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(SAMPLES.index), # nur samples UA DEM SAMPLES.TSV SIND GÜLTIG
    unit="|".join(UNITS["unit"]), # nur units aus dem units.tsv sind gültig




RDS = ['R1','R2']

REFERENCE = config["ref"]["folder"]
FASTA = config["ref"]["fasta"]
CHLOROPLAST = config["chloroplast"]

CONTEXTS = config["methylation_contexts"]




BISMARK_GENOME_PREPARATION = config["bismark"]["bismark_genome_preparation"]
BISMARK = config["bismark"]["bismark"]
DEDUPLICATE_BISMARK = config["bismark"]["deduplicate_bismark"]
BISMARK_METHYLATION_EXTRACTOR = config["bismark"]["bismark_methylation_extractor"]
BAM2NUC = config["bismark"]["bam2nuc"]
COVERAGE2CYTOSINE = config["bismark"]["coverage2cytosine"]
BISMARK2SUMMARY =  config["bismark"]["bismark2summary"]


VIEWBS = config["viewbs"]["viewbs"]

CREATEVIEWBSFILELIST = config["scripts"]["createviewbsfilelist"]
BISMARKCOV2TILES = config["scripts"]["bismarkcov2tiles"]


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = UNITS.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_trimmed_fastq(wildcards):
	"""Get trimmed files of given sample."""
	fin1 = []
	fin2 = []
	unit = UNITS[UNITS["sample"] == wildcards.sample]["unit"].values
	for i in unit:
		fin1.append("1_trimmed_reads/{}-{}.1.fastq.gz".format(wildcards.sample,i))
		fin2.append("1_trimmed_reads/{}-{}.2.fastq.gz".format(wildcards.sample,i))
	return {"fin1": fin1, "fin2": fin2}

def get_bam(wildcards):
	"""Get bam files for all units of given sample."""
	fin1 = []
	unit = UNITS[UNITS["sample"] == wildcards.sample]["unit"].values
	for i in unit:
		fin1.append("2_bismark_mapped/{}-{}.1_bismark_bt2_pe.bam".format(wildcards.sample,i))
	return {"bam": fin1}

def get_trimmed_fastq1(wildcards):
    """Get trimmed files of given sample."""
    fin1 = []
    unit = UNITS[UNITS["sample"] == wildcards.sample]["unit"].values
    for i in unit:
        fin1.append("1_trimmed_reads/{}-{}.1.fastq.gz".format(wildcards.sample,i))

    return  fin1

def get_trimmed_fastq2(wildcards):
    """Get trimmed files of given sample."""
    fin2 = []
    unit = UNITS[UNITS["sample"] == wildcards.sample]["unit"].values
    for i in unit:
        fin2.append("1_trimmed_reads/{}-{}.2.fastq.gz".format(wildcards.sample,i))
    return  fin2


# die namen der rules muss es auch geben

#ruleorder: trim_galore > mergeFastqR1 > mergeFastqR2

rule all:
    input:
        expand("4_methylation_extraction/{sample}_{context}_100bp.csv", sample = SAMPLES.index, context=CONTEXTS),
	expand("4_methylation_extraction/{sample}_{context}.contextWithStrand.cov.gz", sample = SAMPLES.index, context=CONTEXTS),
	expand("3_deduplicated/{sample}_pe.deduplicated.bam", sample = SAMPLES.index),
        #expand("1_trimmed_reads/{sample}-{unit}.{rd}_fastqc.html",sample = SAMPLES.index, unit = UNITS["unit"] , rd = [1,2]),
        "6_viewBS/methCoverage",
        "6_viewBS/methGlobal",
        "6_viewBS/MethLevDist",
        expand("6_viewBS/MethGeno_{context}", context = CONTEXTS),
        "6_viewBS/BisNonConvRate",
	"5_multiqc",
	expand("bam2nuc/{sample}_pe.deduplicated.nucleotide_stats.txt", sample = SAMPLES.index)



rule trimmomatic_pe:
    input:
        unpack(get_fastq)
    output:
        r1="1_trimmed_reads/{sample}-{unit}.1.fastq.gz",
        r2="1_trimmed_reads/{sample}-{unit}.2.fastq.gz",
        r1_unpaired=temp("1_trimmed_reads/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("1_trimmed_reads/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog=temp("1_trimmed_reads/{sample}-{unit}.trimlog.txt")
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    wrapper:
        "0.27.1/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input:
        fq = "1_trimmed_reads/{sample}-{unit}.{rd}.fastq.gz"
    output:
        html="1_trimmed_reads/{sample}-{unit}.{rd}_fastqc.html",
        zip="1_trimmed_reads/{sample}-{unit}.{rd}_fastqc.zip"
    log:
        "logs/fastqc_trimmed_{sample}_{unit}_{rd}.log"
    message:
        "---------- running fastqc for trimmed fastq files  ----------"
    shell:
        "fastqc {input.fq}"



rule bismark_genome_preparation:
    input:
        REFERENCE
    output:
        REFERENCE+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        REFERENCE+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        #pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 ",
        verbose = "--verbose "
    log:
        'logs/bismark_genome_preparation.log'
    message: """ --------  converting Genome into Bisulfite analogue ------- """
    shell:
        "{BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"


rule samtools_index:
    input: FASTA
    output: FASTA+".fai"
    log:
        'logs/samtools_index.log'
    message:
        """ ---------- indexing reference file ------------ """
    shell:
        "samtools faidx {input}"





rule bismark_alignment:
    input:
        fin1 = "1_trimmed_reads/{sample}-{unit}.1.fastq.gz",
        fin2 = "1_trimmed_reads/{sample}-{unit}.2.fastq.gz",
        CT_conv = REFERENCE+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        GA_conv = REFERENCE+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    output:
        "2_bismark_mapped/{sample}-{unit}.1_bismark_bt2_pe.bam",
        "2_bismark_mapped/{sample}-{unit}.1.fastq.gz_unmapped_reads_1.fq.gz",
        "2_bismark_mapped/{sample}-{unit}.2.fastq.gz_unmapped_reads_2.fq.gz",
#        "2_bismark_mapped/{sample}_pe.nucleotide_stats.txt",
        "2_bismark_mapped/{sample}-{unit}.1_bismark_bt2_PE_report.txt"
    threads: 4
    params:
#        N = "-N 1",
#        L = "-L 20",
        X = "-X 600",
        genomeFolder = "--genome_folder " + REFERENCE,
        outdir = "--output_dir  "+ "2_bismark_mapped/",
        tempdir = "--temp_dir " + "2_bismark_mapped/",
        unmapped = "--un",
        useBowtie2  = "--bowtie2 ",
    log:
        "logs/{sample}_{unit}_bismark_pe_mapping.log"
    message: """-------------   Mapping reads to genome for {wildcards.sample} {input.fin1} {input.fin2}. ------------- """
    shell:
    	"{BISMARK} {params} --multicore {threads}  -1 {input.fin1} -2 {input.fin2} 2> {log}"




rule deduplicate_bismark:
    input:
        unpack(get_bam)
    output:
        files = "3_deduplicated/{sample}_pe.deduplicated.bam"
    params:
        pe = "-p",
        outdir = "-o 3_deduplicated",
        bam = "--bam",
        mult = "--multiple",
        mv = "3_deduplicated/{sample}-*.bam"
    log:
        "logs/deduplicate_bismark_{sample}.log"
    message: """-----------    deduplicating reads of {input.bam} to {output.files} ----------"""
    shell:
        "{DEDUPLICATE_BISMARK} {params.pe} {params.outdir} {params.bam} {params.mult}  {input.bam} && mv {params.mv} {output.files}"



rule bam2nuc:
    input:
        fa=REFERENCE,
        bam="3_deduplicated/{sample}_pe.deduplicated.bam"
    output:
        "bam2nuc/{sample}_pe.deduplicated.nucleotide_stats.txt"
    log:
        "logs/bam2nuc_{sample}.log"
    message: """ -------- calculating nucleotide stats for {input.bam}  ------- """

    shell:
        "{BAM2NUC} --dir bam2nuc --genome_folder {input.fa}  {input.bam} "

rule bismark_methylation_extractor:
    input:
        fa=REFERENCE,
        file="3_deduplicated/{sample}_pe.deduplicated.bam"
    output:
        "4_methylation_extraction/CpG_context_{sample}_pe.deduplicated.txt.gz",
        "4_methylation_extraction/CHG_context_{sample}_pe.deduplicated.txt.gz",
        "4_methylation_extraction/CHH_context_{sample}_pe.deduplicated.txt.gz",
        "4_methylation_extraction/{sample}_pe.deduplicated_splitting_report.txt",
        "4_methylation_extraction/{sample}_pe.deduplicated.M-bias.txt",
        "4_methylation_extraction/{sample}_pe.deduplicated.bedGraph.gz",
        "4_methylation_extraction/{sample}_pe.deduplicated.bismark.cov.gz"
#	"4_methylation_extraction/{sample}.bis_rep.cov.CX_report.txt"
    threads: 30
    params:
        outdir = " 4_methylation_extraction "
    log:
        "logs/{sample}.bismark_methylation_extractor.log"
    message: """ -------- extracting methylation information for {input.file}  ------- """

    shell:
        "{BISMARK_METHYLATION_EXTRACTOR} --comprehensive --ignore_r2 2 --gzip --multicore 10 --bedGraph --CX --cytosine_report --genome_folder {input.fa} -o {params.outdir} --buffer_size 10G {input.file}   2> {log}"


rule genomewidecytosinemethylationreport:
    input:
        fa=REFERENCE,
        bismarkcov="4_methylation_extraction/{sample}_pe.deduplicated.bismark.cov.gz"
    output:
        "4_methylation_extraction/{sample}.bis_rep.cov.CX_report.txt"
    log:
     "logs/bismark_genomewidecytosinemethylationreport_{sample}.log"
    message: """ -------- Witing genome wide cytosine methylation report for {input.bismarkcov}  ------- """
    shell:
        "{COVERAGE2CYTOSINE} -CX -o {output} --genome_folder {input.fa} {input.bismarkcov}"

rule summarizedinucleotides:
    input:
        fa=REFERENCE,
        bismarkcov="4_methylation_extraction/{sample}_pe.deduplicated.bismark.cov.gz"
    output:
        "4_methylation_extraction/{sample}.bis_rep.cov.DN_report.txt.CpG_report.txt",
        "4_methylation_extraction/{sample}.bis_rep.cov.DN_report.txt.CpG_report.merged_CpG_evidence.cov"
    log:
        "logs/bismark_summarizedinucleotides_{sample}.log"
    message: """ -------- Witing genome wide cytosine methylation report with merged CpG for {input.bismarkcov}  ------- """
    shell:
        "{COVERAGE2CYTOSINE} --merge_CpG -o {output} --genome_folder {input.fa} {input.bismarkcov}"



#rule cContextCovFileFromGWCMR:
#    input:
#        gwcmr="4_methylation_extraction/{sample}.bis_rep.cov.CX_report.txt"
#    output:
#        "4_methylation_extraction/{sample}_{context}.context.cov.gz"
#    threads: 3
#    params:
#        con = "{context}"
#    log:
#        "logs/MethylationCoverageFiles/{sample}_{context}.context.cov.log"
#    message: """ -------- Writing methylation coverage file for sample {wildcards.sample} and context {wildcards.context}  ------- """
#    shell:
#        " grep "r'"\s{params.con}\s"'" {input.gwcmr}   "
#        "| awk   '{{ if ($4 == 0 && $5 == 0) ; else if ($4 == 0) "
#        "  printf ("r'"%s\t%s\t%s\t%s\t%s\t%s\n"'", $1, $2, $2, 0, $4, $5);"
#        " else  printf ("r'"%s\t%s\t%s\t%s\t%s\t%s\n"'", $1, $2, $2, ($4/($4+$5))*100, $4, $5);}}' "
#        "| gzip > {output}"


rule cContextCovFileWithStrandInfoFromGWCMR:
    input:
        gwcmr="4_methylation_extraction/{sample}.bis_rep.cov.CX_report.txt"
    output:
        "4_methylation_extraction/{sample}_{context}.contextWithStrand.cov.gz"
    threads: 3
    params:
        con = "{context}"
    log:
        "logs/MethylationCoverageFiles/{sample}_{context}.contextWithStrand.cov.log"
    message: """ -------- Writing methylation coverage file with strand info for sample {wildcards.sample} and context {wildcards.context}  ------- """
    shell:
        " grep "r'"\s{params.con}\s"'" {input.gwcmr}   "
        "| awk   '{{ if ($4 == 0 && $5 == 0) printf ("r'"%s\t%s\t%s\tNA\t%s\t%s\t%s\n"'", $1, $2, $2, $4, $5, $3) ; else if ($4 == 0) "
        "  printf ("r'"%s\t%s\t%s\t%s\t%s\t%s\t%s\n"'", $1, $2, $2, 0, $4, $5, $3);"
        " else  printf ("r'"%s\t%s\t%s\t%s\t%s\t%s\t%s\n"'", $1, $2, $2, ($4/($4+$5))*100, $4, $5, $3);}}' "
        "| gzip > {output}"
												    


rule cytosinecoveragetiling: # einen pro C context
    input:
        "4_methylation_extraction/{sample}_{context}.contextWithStrand.cov.gz"
    output:
        "4_methylation_extraction/{sample}_{context}_100bp.csv"
    log:
        "logs/cytosinecoveragetiling_{sample}_{context}.log"
    message: """ -------- merging coverage for sample {wildcards.sample} and context {wildcards.context} in tiles   ------- """
    shell:
        "zcat {input} | python {BISMARKCOV2TILES} > {output}"


rule bgzipgwcmr:
    input:
        "4_methylation_extraction/{sample}.bis_rep.cov.CX_report.txt"
    output:
        "6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz"
    log:
        "logs/bgzipgwcmr_{sample}.log"
    message: """ -------- bgzip genome wide methylation coverage file for {input}  ------- """
    shell:
        "bgzip -c  {input} > {output}"


rule csigwcmr:
    input:
        "6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz"
    output:
        "6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi"
    log:
        "logs/csigwcmr_{sample}.log"
    message: """ -------- tabix csi indexing genome wide methylation coverage file for {input}  ------- """
    shell:
        "tabix -p vcf --csi {input}"


rule ViewBS_samplefile:
    input:
        files=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz", sample=SAMPLES.index)
    output:
        "6_viewBS/inputfiles.txt"
    log:
        "logs/ViewBS_samplefile.log"
    message: """ -------- create list of samples for ViewBS  ------- """
    shell:
        "python {CREATEVIEWBSFILELIST} {input.files} > 6_viewBS/inputfiles.txt"

rule ViewBS_MethCoverage:
    input:
        files="6_viewBS/inputfiles.txt",
        fa=FASTA,
	fai = FASTA + ".fai",
        indexes=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi", sample=SAMPLES.index)

    output:
        directory("6_viewBS/methCoverage")
    log:
        "logs/ViewBS_MethCoverage.log"
    message: """ -------- run ViewBS meth coverage ------- """
    shell:
        "{VIEWBS} MethCoverage --reference {input.fa} --sample file:{input.files}  --outdir {output}  "


rule ViewBS_GlobalMethLev:
    input:
        files="6_viewBS/inputfiles.txt",
        indexes=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi", sample=SAMPLES.index)

    output:
        directory("6_viewBS/methGlobal")
    log:
        "logs/ViewBS_GlobalMethLev.log"
    message: """ -------- run ViewBS GlobalMethLev ------- """
    shell:
        "{VIEWBS} GlobalMethLev --sample file:{input.files}  --outdir {output}  "


rule ViewBS_MethLevDist:
    input:
        files="6_viewBS/inputfiles.txt",
        indexes=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi", sample=SAMPLES.index)

    output:
        directory("6_viewBS/MethLevDist")
    params:
        bin = " --binMethLev 0.1"
    log:
        "logs/ViewBS_MethLevDist.log"
    message: """ -------- run ViewBS methLevDist ------- """
    shell:
        "{VIEWBS} MethLevDist {params} --sample file:{input.files}  --outdir {output}  "


rule ViewBS_MethGeno:
    input:
        files="6_viewBS/inputfiles.txt",
        fa=FASTA,
	    fai = FASTA + ".fai",
        indexes=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi", sample=SAMPLES.index)
    output:
        directory("6_viewBS/MethGeno_{context}")
    params:
        genomeLength = "--genomeLength " + FASTA + ".fai",
        context = "--context {context}"
    log:
        "logs/ViewBS_MethGeno_{context}.log"
    message: """ -------- run ViewBS MethGeno ------- """
    shell:
        "{VIEWBS} MethGeno {params} --sample file:{input.files}  --outdir {output}  "

rule ViewBS_BisNonConvRate:
    input:
        files="6_viewBS/inputfiles.txt",
        indexes=expand("6_viewBS/{sample}.bis_rep.cov.CX_report.txt.gz.csi", sample=SAMPLES.index)
    output:
        directory("6_viewBS/BisNonConvRate")
    params:
        chloroplast = "--chrom " + CHLOROPLAST
    log:
        "logs/BisNonConvRate.log"
    message: """    ----------  run ViewBS BisNonConvRate   ----------   """
    shell:
        "{VIEWBS} BisNonConvRate {params} --sample file:{input.files} --outdir {output} "



#rule bismark_reports:
#    input:     
#        expand(["3_deduplicated/{sample}-1.1_bismark_bt2_pe.multiple.deduplication_report.txt","4_methylation_extraction/{sample}_pe.deduplicated.M-bias.txt","4_methylation_extraction/{sample}_pe.deduplicated_splitting_report.txt"],sample = SAMPLES.index)
#    output:
#        expand(["5_bismark_summary/{sample}-1.1_bismark_bt2_pe.multiple.deduplication_report.txt", "5_bismark_summary/{sample}_pe.deduplicated.M-bias.txt", "5_bismark_summary/{sample}_pe.deduplicated_splitting_report.txt"],sample=SAMPLES.index)
#    log:
#        "logs/bismark_reports.log"
#    message: """    ----------  bismark_reports  ----------   """
#    shell:
#        "rm -r 5_bismark_summary &&"
#        "mkdir 5_bismark_summary && "
#        "cd 5_bismark_summary && "
#        "ln -s ../2_bismark_mapped/*.bam ../2_bismark_mapped/*.txt ../3_deduplicated/*.txt ../4_methylation_extraction/*M-bias.txt ../4_methylation_extraction/*splitting_report.txt ."
        

#rule bismark_summary_report:
#    input:
#        expand(["5_bismark_summary/{sample}-1.1_bismark_bt2_pe.multiple.deduplication_report.txt", "5_bismark_summary/{sample}_pe.deduplicated.M-bias.txt", "5_bismark_summary/{sample}_pe.deduplicated_splitting_report.txt"],sample=SAMPLES.index)
#    output:
#        "5_bismark_summary/bismark_summary_report.txt"
#    log:
#        "logs/bismark_summary_report.log"
#    message: """ -----------  running bismark_summary_report for {input}  -------------"""
#    shell:
#        "cd 5_bismark_summary && {BISMARK2SUMMARY} "


rule multiqcSummaryReport:
    input:
        expand(["3_deduplicated/{sample}-1.1_bismark_bt2_pe.multiple.deduplication_report.txt","4_methylation_extraction/{sample}_pe.deduplicated.M-bias.txt","4_methylation_extraction/{sample}_pe.deduplicated_splitting_report.txt"],sample = SAMPLES.index)
    output:
        directory("5_multiqc")
    log:
        "logs/multiqc.log"
    message: """ -----------  running multiqc  -------------"""
    shell:
        "multiqc -o {output} . "
