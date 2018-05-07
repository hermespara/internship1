import subprocess

########################################################
############								############
############	   SRAtookit parameters		############
############                                ############
########################################################
# Indicate workdir. Normally the ncbi directory.
# This directory has to contain the Snakefile
WorkDir = '/data/home/hparaqindes/ncbi/'
# dbGaP : the project name where sra data will be stored.
# The directory should have the same name as the project. 
dbGaP = 'dbGaP-13871'
# public : by default sratoolkit directory.
# In this directory will be stored all refseq for the sra files
public = 'public'
# sra_toolkit_Dir : the path to the sratoolkit directory. 
# Notice the path should end to bin directory
# I suggest to use the 2.8.1-1 version as i had some problems with prefetch and aspera with the newest version
# See the readme for more information
sra_toolkit_Dir = '/data/home/hparaqindes/sratoolkit.2.8.1-1-ubuntu64/bin/'
# vdb_validate, fastq_dump and prefetch don't need to be change. 
# Exept in case of un update of sratoolkit
vdb_validate = 'vdb-validate'
fastq_dump = 'fastq-dump'
prefetch = 'prefetch'
# aspera : path to executable aspera and asperaweb_id_dsa.openssh
# Example "<path/to/.aspera/connect/bin/ascp|<path/to/.aspera/connect/etc/asperaweb_id_dsa.openssh"
aspera = "\"/data/home/hparaqindes/.aspera/connect/bin/ascp|/data/home/hparaqindes/.aspera/connect/etc/asperaweb_id_dsa.openssh\""
# SAMPLES variable is a list. It will contain all the sra run numbers.
SAMPLES = []
# sra_id :  path to the file that contains all the sra run numbers listed in a column.
sra_id = '/data/home/hparaqindes/ncbi/SRA_numbers_try'

# open the sra_id and append to the SAMPLES list all the sra run ID.
with open(sra_id, "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
# print the list of sra run ID
print(SAMPLES)


########################################################
############								############
############	    FASTQC parameters		############
############                                ############
########################################################
fastqc_software = '/data/home/hparaqindes/FastQC/fastqc'

#########################################################################
######														#############
######				STAR parameters							#############
######														#############
#########################################################################

#input the Star directory
star = '/newdisk/tals/STAR-2.5.0b/bin/Linux_x86_64/STAR'
# the genomes Directory used by STAR
genome_Dir = '/newdisk/tals/IRFinder-IRFinder-1.2.0/Human_hg19_GRCh37_75_ref_with_blacklist_intronLength_1_ROI_nonPolyA/STAR'


#########################################################################
######														#############
######				IRF parameters							#############
######														#############
#########################################################################

IRF_genome = '/newdisk/tals/IRFinder-IRFinder-1.2.0/Human_hg19_GRCh37_75_ref_with_blacklist_intronLength_1_ROI_nonPolyA'
IRF_dir = '/newdisk/tals/IRFinder-IRFinder-1.2.0/bin/IRFinder'

########################################################
############								############
############	    JunctionCover2IRF		############
############                                ############
########################################################
junctionCover2IRF_file_dir = '/data/home/hparaqindes/junctionsCover2IRF.py'
introns_file = '/newdisk/tals/IRFinder/IRFinder.U12_ID_geneName_pos_strand.not_analyzed.txt'

########################################################
############								############
############	    BAM anonymize			############
############                                ############
########################################################
anonymize_bam_file = '/data/home/hparaqindes/anonymize_bam_file.py'

########################################################
############								############
############	    HTSeq Count 			############
############                                ############
########################################################

htseq_count_gtf = '/newdisk/tals/Annotations/Homo_sapiens.GRCh37.75.gtf'
# I installed the new version of htseq-count in my home.
# pip install --user HTSeq
# or upgrade pip install --user HTSeq --upgrade if you have an older version
htseq_count_software = '/data/home/hparaqindes/.local/bin/htseq-count'


########################################################
############								############
############	    Snakemake rules			############
############                                ############
########################################################

########################################################
############								############
############	        Rule all			############
############                                ############
########################################################

# Checks the final files. If one of the files is missing,
# snakemake will show an error. Each one of this files,
# should be an output of the other rules. Notice that 
# some files are commented "#" because they are deleted
# during the pipeline. If you want to keep these files,
# uncomment the line and modify the rule where they are
# deleted.

rule all:
    input:
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.zip', sample = SAMPLES),
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam', sample=SAMPLES),
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + 'IRFinder-IR-nondir.txt', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + '{sample}', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam.bai', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}.count', sample = SAMPLES)

########################################################
############								############
############     Rule dowload_sra_files		############
############                                ############
########################################################

# This rule will download all sra files using prefetch of sratoolkit and aspera.
# It will move the files into a directory with the same name as the basename.
# Exemple SRR12345.sra will be moved to the directory SRR12345
# You need to have installed aspera.
# It takes time to download. 
# The prefetch command needs to be executed in the dbGaP directory.

rule download_sra_files:
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra'
    shell:
        """
        cd {WorkDir}{dbGaP}
        {sra_toolkit_Dir}{prefetch} -a {aspera} {wildcards.sample}
        mv {WorkDir}{dbGaP}/sra/{wildcards.sample}.* {WorkDir}{dbGaP}/sra/{wildcards.sample}/
        """
########################################################
############								############
############	    Rule fastq_dump			############
############                                ############
########################################################

# This rule will use fastq-dump of sratoolkit to tranform the sra file to fastq.gz
# There are two output files: {sample}_1.fastq.gz and {sample}_2.fastq.gz.
# Notice that in this rule, the sra file is deleted at the end. If you want
# to keep this file, please remove or echo the last line of shell command.
# –gzip Compress output using gzip
# –split-3 will output 1,2, or 3 files: 1 file means the data is not paired. 2 files means paired data with no low quality reads or reads shorter than 20bp. 3 files means paired data, but asymmetric quality or trimming.

rule fastq_dump:
	input:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra'
	output:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
	shell:
		"""
        cd {WorkDir}{dbGaP}/sra/{wildcards.sample}
		pwd
		{sra_toolkit_Dir}{fastq_dump} --gzip --split-3 {input}
        echo removing {wildcards.sample}.sra
        rm {wildcards.sample}.*
		"""

########################################################
############								############
############	    Rule fastQC				############
############                                ############
########################################################

# This rule uses fastqc software to check the quality of fastq files. 
# In input it takes the fastq files generated from fastq-dump.
# it outputs 4 files in a directory called fastQC.
# This rule can be runed before or after the star_mapping rule. 
# -t number of threads used
# -o output direcotry.

rule fastQC:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.html',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.zip',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.html',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.zip'
    shell:
        """
        {fastqc_software} -t 8 -o {WorkDir}{dbGaP}/sra/{wildcards.sample}/fastQC {input[0]} {input[1]}
        """
########################################################
############								############
############	    Rule star_mapping		############
############                                ############
########################################################

# This rule will use STAR software to mapp the files to the human index
# It takes the fastq files as input and all the output will be stored
# in a directory called STAR. All the files will have as prefix the sra
# number. The bam files will be deleted at the bam_anonymize rule. If 
# you want to keep the bam files, please modify bam_anonymize rule.
# --runThreadN Number of threads
# --limitBAMsortRAM the max of RAM allowed to use
# --readFilesCommand : UncompressionCommand. Here we have .gz files so the zcat option is added.
# --genomeDir the genome directory
# --outFilterMultimapNmax : max  number  of  multiple  alignments  allowed  for  a  read.
# If  exceeded,  the  read  is  considered unmapped. Here 1.
# --outSAMunmapped : output of unmapped reads in the SAM format. By default: None. Here None
# --outSAMtype : type of SAM/BAM output. Here BAM and Unsorted for 2 argument and SortedByCoordinate for the 3.
# --readFilesIn : input files. Here the fastq files
# --outFileNamePrefix : The output files prefix. Here the sra number
 
rule star_mapping:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.out.bam'
    shell:
        """
        {star} --runThreadN 8 \
        --limitBAMsortRAM 30000000000 \
        --readFilesCommand zcat \
        --genomeDir {genome_Dir} \
        --outFilterMultimapNmax 1 \
        --outSAMunmapped None \
        --outSAMtype BAM Unsorted \
        SortedByCoordinate \
        --readFilesIn {input[0]} {input[1]} \
        --outFileNamePrefix {WorkDir}{dbGaP}/sra/{wildcards.sample}/STAR/{wildcards.sample}_
		"""

########################################################
############								############
############	    Rule samtools_index		############
############                                ############
########################################################

# here we use the samtools index software. To be executed, it needs the bam file as input.
# Notice that in input there are fastq files also. These two files will be deleted at the end.
# If you want to keep it, modify the shell command. 
# The index will output a .bam.bai file in the STAR directory. 
# The index is nedded for the bam_anonymize rule as pysam needs the .bai file

rule samtools_index:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input[0]}
        echo removing {input[1]}
        rm {input[1]}
        echo removing {input[2]}
        rm {input[2]}
        """

########################################################
############								############
############	    Rule run_IRF			############
############                                ############
########################################################

# The run_IRF rule uses the IRFinder software to detect intron retention.
# It takes as input the {sample}_Aligned.sortedByCoord.out.bam
# All the output files will be stocked in IRFinder directory
# -m input file format
# -r IRF genome direcotory
# -d output directory

rule run_IRF:
	input:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam'
	output:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + 'IRFinder-IR-nondir.txt'
	shell:
		"""
		{IRF_dir} -m BAM -r {IRF_genome} \
		-d {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder {input}
		"""

########################################################
############								############
############	    Rule run_IRF			############
############                                ############
########################################################

# This rule get the junctions counts for all specify introns and add them
# to an IRFinder file.
# It takes the {sample}_Aligned.sortedByCoord.out.bam as input.
# Output a file with the same name as the sra number in IRFinder directory
# Uses python2.7
# --OutFolder : output Folder
# --AlignmentPrefix : the prefix of STAR files
# --Introns : file that contains all the introns
# --TMP : temporary directory.
# Notice that the junctionCover2IRF.py uses samtools view
# The output:
# 1chr	2:start(0-based)	3:end	4:GeneName/GeneID/staticWarning	5:0	6:Strand	7:ExclBases(=intron length)	8:Coverage(0)
# 9:IntronDepth(0)	10:Intron25thPercentile(0)	11:Intron50thPercentile(0) 12:Intron75thPercentile(0)	
# 13:ExonIntronJunctionLeft(samtools view)	14:ExonIntronJunctionRight(samtools view)	15:IntronDepthFirst50(0)	16:IntronDepthLast50(0)	
# 17:SpliceLeft(Splice, samtools view)	18:SpliceRight(Splice, samtools view)	19:SpliceExact(samtools view)	20:IRratio(0)	21:Warnings


rule junctionCover2IRF:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + '{sample}'
    shell:
        """
        python2.7 {junctionCover2IRF_file_dir} \
        --OutFolder {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder \
		--AlignmentPrefix {WorkDir}{dbGaP}/sra/{wildcards.sample}/STAR/{wildcards.sample} \
        --Introns {introns_file} \
        --startCorection \
		--TMP {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder
        """

########################################################
############								############
############	    Rule htseq_count		############
############                                ############
########################################################

# This rule uses htseq-count software. It counts for each gene how many reads map to it.
# As input it takes {sample}_Aligned.sortedByCoord.out.bam file and the gtf file
# It outputs an {sample}.count file at STAR directory
# -f : input format file.
# -r : Order. For paired-end data, the alignment have to be sorted either by read name or by alignment position. Here pos
# --nonunique none : --nonunique none (default): the read (or read pair) is counted as ambiguous and not counted for any features. Also, if the read (or read pair) aligns to more than one location in the reference, it is scored as alignment_not_unique.
# -m : Mode to handle reads overlapping more than one feature. Possible values for <mode> are union, intersection-strict and intersection-nonempty
# --additional-attr : Additional feature attributes, which will be printed as an additional column after the primary attribute column but before the counts column(s). The default is none, a suitable value to get gene names using an Ensembl GTF file is gene_name.
# -t : feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon)
# -i : GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is gene_id.

rule htseq_count:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}.count'
    shell:
        """
        {htseq_count_software} -f bam \
        -r pos \
        -s no \
        --nonunique none \
        -m intersection-nonempty \
        --additional-attr gene_name \
        -t exon \
        -i gene_id \
        {input} {htseq_count_gtf} > {output}
        """

########################################################
############								############
############	    Rule bam_anonymize		############
############                                ############
########################################################

# This rule will anonymize the bam file outputed by star.
# It uses python 2.7 with pysam library.
# Input : {sample}_Aligned.sortedByCoord.out.bam
# Output : {sample}_Aligned_anonymize.sortedByCoord.out.bam in the STAR_anonymized directory
# Notice that this rule will also delete the bam and bai file in the STAR directory. 
# If you want to keep the files, modify the shell command.
# The other input files, indicate to snakemake, that this rule needs to be run after the htseq_count rule is finished;
# --bamIn : input bam file
# --bamOut :  outut bam file

rule bam_anonymize:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.out.bam',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}.count'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam'
    shell:
        """
        python2.7 {anonymize_bam_file} \
        --bamIn {input[0]} \
        --bamOut {output}
        echo removing {input[0]}
        rm {input[0]}
        echo removing {input[1]}
        rm {input[1]}
        echo removing {input[2]}
        rm {input[2]}
        """

########################################################
############								############
############ Rule new_bam_anonymized_index	############
############                                ############
########################################################

# This rule will index the anonymized bam file.
# Input : {sample}_Aligned_anonymize.sortedByCoord.out.bam
# Output : {sample}_Aligned_anonymize.sortedByCoord.out.bam.bai
# It can be used to be visualized in IGV

rule new_bam_anonymized_index:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input}
        """

