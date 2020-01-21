fastq_files, = glob_wildcards('data/02_seqdata/{filename}.fq.gz')
  
rule all:
  input: 
    expand("results/fastqc/{filename}_fastqc.zip", filename = fastq_files)

rule fastqc:
  input: 
    "data/02_seqdata/{sample}.fq.gz"
  output:
    html="results/fastqc/{sample}.html",
    zip="results/fastqc/{sample}_fastqc.zip"
  log:
    "logs/fastqc/{sample}.log"
  threads: 
    1
  shell:
    '''
    nice --adjustment=+30 fastqc {input} -o=results/fastqc/ 2> {log}
    '''
