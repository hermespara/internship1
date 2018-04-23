import subprocess


WorkDir = '/home/hermesparaqindes/Bureau/ncbi/'
dbGaP = 'dbGaP-13871'
public = 'puplic'
sra_toolkit_Dir = '/home/hermesparaqindes/Téléchargements/sratoolkit.current-ubuntu64/sratoolkit.2.9.0-ubuntu64/bin/'
vdb_validate = 'vdb-validate'
fastq_dump = 'fastq-dump'
prefetch = 'prefetch'
aspera = "\"/home/hermesparaqindes/.aspera/connect/bin/ascp|/home/hermesparaqindes/.aspera/connect/etc/asperaweb_id_dsa.openssh\""
SAMPLES = []
dr = '/home/hermesparaqindes/Bureau/Snake_make/'
sra_id = '/home/hermesparaqindes/Bureau/Snake_make/SRA_numbers_try'
#sra = expand('{sra}', sra=sra_id)
with open('/home/hermesparaqindes/Bureau/Snake_make/SRA_numbers_try', "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
print(SAMPLES)
rule all:
    input:
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES),
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.html', sample = SAMPLES)
        #WorkDir+'{sample}._1.gz'
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
        expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}.sra', sample = SAMPLES)
    shell:
        """
        cd {WorkDir}
        ./bashFile.pbs {input.r1}
        """

rule fastq_dump_:
    input:
        outdir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}',sample = SAMPLES),
        file_in = expand(WorkDir + dbGaP + "/sra" + "/" +'{sample}'+'/{sample}.sra', sample = SAMPLES)

    output:
        #directory = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}',sample = SAMPLES),
        out1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES),
        out2 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_2.fastq.gz', sample = SAMPLES)
    shell:
        """
        cd {input.outdir}
        fastq-dump --gzip --split-3 {input.file_in}
        """
rule fastQC:
    input:
        fastqDir = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}', sample = SAMPLES),
        fastq_1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1.fastq.gz', sample = SAMPLES)
    output:
        html_1 = expand(WorkDir + dbGaP + "/sra" + "/" + '{sample}'+'/{sample}_1_fastqc.html', sample = SAMPLES)
    shell:
        """
        cd {input.fastqDir}
        fastqc -t 8 {input.fastq_1}
        """

