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
STAR = '/home/hermesparaqindes/ncbi-outdir/STAR-2.5.2a/bin/Linux_x86_64/STAR'
genome_Dir = '/home/hermesparaqindes/....'
#sra = expand('{sra}', sra=sra_id)
with open('/home/hermesparaqindes/Bureau/Snake_make/SRA_numbers_try', "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
print(SAMPLES)
a = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES)
print(a)

rule all:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra.vdbcache.cache', sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + '/{sample}.sra.vdbcache.cache', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.html', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.zip', sample=SAMPLES)
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
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra.vdbcache.cache', sample = SAMPLES)
    shell:
        """
        cd {WorkDir}
        ./bashFile.pbs {input.r1}
        """

rule fastq_dump_:
    input:
        outdir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample=SAMPLES),
        file_in = expand(WorkDir + dbGaP + "/sra" + "/" +'{sample}'+'/{sample}.sra', sample=SAMPLES)
    output:
        #directory = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}',sample = SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + '/{sample}.sra.vdbcache.cache', sample=SAMPLES)
    shell:
        """
        cd {input.outdir}
        pwd
        echo {input.file_in}
        fastq-dump --gzip --split-3 {input.file_in}
        """

rule fastQC:
    input:
        fastqDir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample = SAMPLES),
        fastq_1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES)
    output:
        html_1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.html', sample = SAMPLES),
        fastqc_zip = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.zip', sample = SAMPLES)
    shell:
        """
        cd {input.fastqDir}
        fastqc -t 8 {input.fastq_1}
        """
'''
rule star_mapping:
    input:
        fastq_gz_1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES),
        fastq_gz_2 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample = SAMPLES)
    output:
        sorted_bal_out = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}' + '')
    shell:
        """
        {STAR} --runThreadN 8 \
        --limitBAMsortRAM 30000000000 \
        --readFilesCommand \
        zcat --genomeDir {genome_Dir} \
        --outFilterMultimapNmax 1 \
        --outSAMunmapped None \
        --outSAMtype BAM Unsorted \
        SortedByCoordinate \
        --readFilesIn {input.fastq_gz_1} {input.fastq_gz_2} \
        --outFileNamePrefix {????}/${????}_
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

