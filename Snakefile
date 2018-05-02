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
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_1_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.html', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/fastQC" + '/{sample}_2_fastqc.zip', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + 'IRFinder-IR-nondir.txt', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + '{sample}', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam.bai', sample = SAMPLES)


rule try_another_way_to_download:
    input:
        r1 = WorkDir + 'SRA_numbers_try',
        r2 = WorkDir + 'bashFile.pbs'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra'
    shell:
        """
        cd {WorkDir}{dbGaP}
        {sra_toolkit_Dir}{prefetch} -a {aspera} {wildcards.sample}
        mv {WorkDir}{dbGaP}/sra/{wildcards.sample}.* {WorkDir}{dbGaP}/sra/{wildcards.sample}/
        """
rule fastq_dump:
	input:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra'
        #expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample=SAMPLES)
	output:
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
		WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
	shell:
		"""
        cd {WorkDir}{dbGaP}/sra/{wildcards.sample}
		pwd
		fastq-dump --gzip --split-3 {input}
		"""

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
        fastqc -t 8 -o {WorkDir}{dbGaP}/sra/{wildcards.sample}/fastQC {input[0]} {input[1]}
        """

rule star_mapping:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam'
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

rule samtools_index:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + "/STAR" + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input}
        """

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


rule junctionCover2IRF:
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + '{sample}'
    shell:
        """
        python2.7 {juncitonCover2IRF_file_dir} \
        --OutFolder {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder \
		--AlignmentPrefix {WorkDir}{dbGaP}/sra/{wildcards.sample}/STAR/{wildcards.sample} \
        --Introns {introns_file} \
        --startCorection \
		--TMP {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder
        """

rule bam_anonymize:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam',
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR' + "/" + '{sample}_Aligned.sortedByCoord.out.bam.bai'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam'
    shell:
        """
        PYTHONPATH=$PYTHONPATH:/usr/remote/python_packages/lib/python2.7/site-packages/
        python2.7 {anonymize_bam_file} \
        --bamIn {input[0]} \
        --bamOut {output}
        """

rule new_bam_anonymized_index:
    input:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam'
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/STAR_anonymized' + "/" + '{sample}_Aligned_anonymize.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input}
        """


'''
rule junctionCover2IRF:
    output:
        WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '/IRFinder' + "/" + '{sample}_Aligned.out.bam'
    shell:
        """
        python2.7 {juncitonCover2IRF_file_dir} \
        --OutFolder {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder \
		--AlignmentPrefix {WorkDir}{dbGaP}/sra/{wildcards.sample}/STAR/{wildcards.sample} \
        --Introns {introns_file} \
        --startCorection \
		--TMP {WorkDir}{dbGaP}/sra/{wildcards.sample}/IRFinder
        """

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

