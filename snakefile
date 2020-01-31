fastq_files, = glob_wildcards('data/02_seqdata/{filename}.fq.gz')

import re

rule all:
  input: 
    multiqc_report = "results/multiqc/multiqc.html",
    salmon_quant = expand("results/salmon/salmon_quant/{sample_id}", sample_id = [re.sub('_[1-2]$', '', string) for string in fastq_files])

rule md5:
  input:
    "data/02_seqdata/MD5.txt",
    "data/02_seqdata/MD5-1.txt"
  log:
    "results/md5/md5_checks.txt"
  output:
    touch("results/md5_checks_completed.done")
  threads: 20
  shell:
    '''
    cat {input} | awk '{{print $1 " data/02_seqdata/" $2}}' | parallel --jobs {threads} --pipe -N1 md5sum -c >> {log}
    '''

rule fastqc:
  input: 
    "results/md5/md5_checks.txt",
    "data/02_seqdata/{fastq_file}.fq.gz"
  output:
    zip="results/fastqc/{fastq_file}_fastqc.zip"
  log:
    "logs/fastqc/{fastq_file}_fastqc.log"
  threads: 
    1
  shell:
    '''
    nice --adjustment=+30 fastqc {input} -o=results/fastqc/ 2> {log}
    '''
    
rule multiqc:
  input: 
    expand("results/fastqc/{filename}_fastqc.zip", filename = fastq_files)
  output: 
    "results/multiqc/multiqc.html"
  shell: 
    '''
    multiqc results/fastqc/ -n {output}
    '''

rule trim_reads:
  input:
    multiqc = "results/multiqc/multiqc.html",
    fastq = "data/02_seqdata/{fastq_file}.fq.gz"
  output: 
    "results/trim_reads/{fastq_file}.fq"
  threads: 1
  shell:
    '''
    seqtk trimfq -b 11 {input.fastq} > {output}
    '''

# Download mouse reference transcriptome from Gencode v23
rule download_mouse:
  output:
    transcriptome = 'data/Mouse/gencode.vM23.transcripts.fa.gz',
    genome = 'data/Mouse/GRCm38.primary_assembly.genome.fa.gz'
  log:
    "logs/download_mouse.txt"
  shell:
    '''
    curl -o {output.transcriptome} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz 2> {log}
    curl -o {output.genome} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz 2>> {log}
    '''

# Index mouse genome for salmon quantification
rule salmon_index:
  input:
    transcriptome = 'data/Mouse/gencode.vM23.transcripts.fa.gz',
    genome = 'data/Mouse/GRCm38.primary_assembly.genome.fa.gz'
  output:
      decoys = temp("results/salmon/mouse_transcriptome_index/decoys.txt"),
      gentrome = temp("results/salmon/mouse_transcriptome_index/gentrome.fa"),
      index = directory("results/salmon/mouse_transcriptome_index/ref_idexing"),
  threads: 30
  params:
      kmer = 31
  shell:
      '''
      # Crate decoy file
      echo "Creating Decoy File"
      grep "^>" <(zcat {input.genome}) | cut -d " " -f 1 > {output.decoys} &&
      sed -i -e 's/>//g' {output.decoys} &&
      # Concatenate genome and transcriptome
      echo "Concatenating genome and transcriptome"
      zcat {input.transcriptome} {input.genome} > {output.gentrome} &&
      # Create index
      echo "Creating index"
      salmon index \
                -t {output.gentrome} \
                -i {output.index} \
                -d {output.decoys} \
                -p {threads} \
                -k {params.kmer} \
                --gencode
      '''

# Quantify read counts via Samon
rule salmon_quant:
  input:
    index = "results/salmon/mouse_transcriptome_index/ref_idexing",
    reads_1 = "results/trim_reads/{sample_id}_1.fq",
    reads_2 = "results/trim_reads/{sample_id}_2.fq"
  output:
    outdir = directory("results/salmon/salmon_quant/{sample_id}")
  params:
    libtype = "ISR",
    numBootstraps = 10,
    minScoreFraction = 0.8,
  threads: 10
  shell:
    '''
    ionice -c2 -n7 \
    salmon quant \
            -i {input.index} \
            -l {params.libtype} \
            -1 {input.reads_1} \
            -2 {input.reads_2} \
            -o {output.outdir} \
            --validateMappings \
            --minScoreFraction {params.minScoreFraction} \
            --numBootstraps {params.numBootstraps}\
            --gcBias \
            --seqBias \
            --writeUnmappedNames \
            --threads {threads} \
    '''
    
