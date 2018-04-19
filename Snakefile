import subprocess


WorkDir = '/home/hermesparaqindes/Bureau/ncbi/sra/'
sra_toolkit_Dir = '/home/hermesparaqindes/Téléchargements/sratoolkit.current-ubuntu64/sratoolkit.2.9.0-ubuntu64/bin/'
vdb_validate = 'vdb-validate'
fastq_dump = 'fastq-dump'
prefetch = 'prefetch'
aspera = '/.aspera/connect/bin/asperaconnect/bin/ascp|/.aspera/connect/bin/asperaconnect/etc/asperaweb_id_dsa.openssh'
SAMPLES = []
dr = '/home/hermesparaqindes/Bureau/Snake_make/'
#sra = expand('{sra}', sra=sra_id)
with open('/home/hermesparaqindes/Bureau/Snake_make/SRA_numbers_try', "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
#print(SAMPLES)
rule all:
    input:
        expand(WorkDir + '{sample}'+'/{sample}.sra', sample=SAMPLES),
        #expand(WorkDir+'{sample}'+'/{sample}_1.fastq.gz', sample=SAMPLES),
        #expand(WorkDir+'{sample}'+'/{sample}_2.fastq.gz', sample=SAMPLES)
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
        r1 = dr + 'SRA_numbers_try',
        r2 = dr + 'bashFile.pbs'
    output:
        WorkDir + '{sample}'+'/{sample}.sra'
    shell:
        """
        cd /home/hermesparaqindes/Bureau/ncbi
        pwd
        echo {input.r2} {input.r1}
        echo {WorkDir}
        {input.r2} {input.r1}
        """

rule fastq_dump_:
    input:
        WorkDir + '{sample}'+'/{sample}.sra'
    output:
        WorkDir + '{sample}'+'/{sample}_1.gz',
        WorkDir + '{sample}'+'/{sample}_2.gz'
    shell:
        """
        cd {WorkDir}{sample}
        {sra_toolkit_Dir}{fastq_dump} --gzip --split-3 {input}
        """
