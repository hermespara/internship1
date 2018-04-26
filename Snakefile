import subprocess


WorkDir = '/home/hermesparaqindes/Bureau/ncbi/'
dbGaP = 'dbGaP-13871'
public = 'puplic'
sra_toolkit_Dir = '/home/hermesparaqindes/Téléchargements/sratoolkit.2.8.1-1-ubuntu64/bin/'
vdb_validate = 'vdb-validate'
fastq_dump = 'fastq-dump'
prefetch = 'prefetch'
aspera = "\"/home/hermesparaqindes/.aspera/connect/bin/ascp|/home/hermesparaqindes/.aspera/connect/etc/asperaweb_id_dsa.openssh\""
SAMPLES = []
dr = '/home/hermesparaqindes/Bureau/Snake_make/'
sra_id = '/home/hermesparaqindes/Bureau/Snake_make/SRA_numbers_try'
star = '/home/hermesparaqindes/ncbi-outdir/STAR-2.5.2a/bin/Linux_x86_64/STAR'
genome_Dir = '/home/hermesparaqindes/....'
#sra = expand('{sra}', sra=sra_id)
with open('/home/hermesparaqindes/Bureau/ncbi/SRA_numbers_try', "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
print(SAMPLES)
a = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES)
print(a)

print(SAMPLES)
#########################################################################
######														#############
######				STAR parameters							#############
######														#############
#########################################################################

#input the Star directory
star = '/home/hermesparaqindes/ncbi-outdir/STAR-2.5.2a/bin/Linux_x86_64/STAR'
# the genomes Directory used by STAR
genome_Dir = '/home/hermesparaqindes/Bureau/Genomes/Human_hg19_GRCh37_75_ref_with_blacklist_intronLength_1_ROI_nonPolyA/STAR'

sample_star_name = expand('{sample}_', sample=SAMPLES)

#########################################################################
######														#############
######				IRF parameters							#############
######														#############
#########################################################################

IRF_output_Dir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/IRFinder", sample = SAMPLES)
IRF_genome = '/home/hermesparaqindes/Bureau/Genomes/Human_hg19_GRCh37_75_ref_with_blacklist_intronLength_1_ROI_nonPolyA'
IRF_dir = '/home/hermesparaqindes/ncbi-outdir/IRFinder-1.2.4/bin/IRFinder'

fastQC_dir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC", sample = SAMPLES)
sample_dir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample = SAMPLES)

rule all:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample = SAMPLES),
		#expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastqc", sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam', sample=SAMPLES)
        #WorkDir+'{sample}._1.gz'
        #/home/hermesparaqindes/Bureau/ncbi/dbGaP-13871/sra/SRR1311916/SRR1311916_1.fastq.gz
'''
rule download_sra_with_prefetch:
    input:
        r1 = dr + '{sample}.txt'
    output:
        dr + '{sample}.sra'
    shell:
        'echo {input.r1} > {output}'
'''
rule try_another_way_to_download:
    input:
        r1 = WorkDir + 'SRA_numbers_try',
        r2 = WorkDir + 'bashFile.pbs'
    output:
    	expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES)
    shell:
        """
        cd {WorkDir}
        ./bashFile.pbs {input.r1}
        """

rule fastq_dump:
	input:
		expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES)
	output:
		expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
		expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES)
	shell:
		"""
		cd {sample_dir}
		pwd
		fastq-dump --gzip --split-3 {input}
		"""

rule fastQC:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample = SAMPLES)
    output:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.zip', sample = SAMPLES)
    shell:
        """
        fastqc -t 8 -o {fastQC_dir} {input[0]} {input[1]}
        """

rule star_mapping:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample = SAMPLES)
    output:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam', sample=SAMPLES)
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
        --outFileNamePrefix {sample_dir}/STAR/{sample_star_name}
		"""
'''
rule samtools_index:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+"/star", sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+"/star"+'{sample}_Aligned.sortedByCoord.out.bam', sample=SAMPLES)
    output:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+"/star", sample=SAMPLES)
    shell:
        """
        cd {input[0]}
        samtools index {input[1]}
        """
rule run_IRF:
    input:

    output:

    shell:
        """
        {IRF_dir} -m BAM -r {IRF_genome} \
        -d {outputDir} {input.sorted_files}
        """
rule junctionCover2IRF:
    input:
        <path/to/junctionsCover2IRF.py>
    output:

'''

