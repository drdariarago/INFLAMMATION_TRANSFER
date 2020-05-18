fastq_files, = glob_wildcards('data/02_seqdata/{filename}.fq.gz')
TISSUE_TYPES = ("placenta_maternal", "lung_maternal", "liver_maternal", "liver_fetal", "placenta_fetal")
MODELS = ("placentas", "fetal_liver", "maternal_liver", "maternal_lung")

MIN_COUNTS = 5
ALPHA = 0.05
MIN_LOGFC = 0.5

import re

rule all:
  input: 
    multiqc_report = "reports/multiqc/multiqc.html",
    linear_model = expand("results/limma_{models}/fitted_model.Rdata", models = MODELS),
    linear_model_summary = expand("results/limma_{models}/fold_change_summary.csv", models = MODELS),
    go_results = expand("results/gost_enrichment_format/enrichment_{models}_{up_or_down}_long.csv",
      models = MODELS, up_or_down = ("upregulated", "downregulated")
    )

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
    "reports/multiqc/multiqc.html"
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
    
# Import sample metadata
rule metadata_import:
  input:
  output:
    rdata = "results/metadata/sample_metadata.Rdata",
    csv = "results/metadata/sample_metadata.csv"
  script:
    "scripts/metadata.R"

# Merge read counts and metadata via tximeta
rule tximeta:
  input:
    salmon_dirs = expand("results/salmon/salmon_quant/{sample_id}", sample_id = list(set([re.sub('_[1-2]$', '', string) for string in fastq_files]))),
    metadata = "results/metadata/sample_metadata.Rdata"
  output:
    expression_results = "results/tximeta/expression_results.Rdata",
    variance_stabilized_counts = "results/tximeta/variance_stabilized_counts.csv"
  script:
    "scripts/tximeta.R"

# Download gene symbols and description
rule download_gene_data:
  input:
    "results/tximeta/variance_stabilized_counts.csv"
  output:
    "results/download_gene_data/gene_names.Rdata"
  script:
    "scripts/download_gene_data.R"

# Fit simple linear models via limma
rule limma:
  input:
    "results/tximeta/expression_results.Rdata"
  output:
    plots = "results/limma/voom_plots.pdf",
    fitted_models = "results/limma/fitted_models.Rdata",
    limma_coefs = "results/limma/limma_coefs.Rdata"
  script:
    "scripts/limma.R"

# Fit linear models for the 2 placentas
rule limma_placentas:
  input:
    expression_results = "results/tximeta/expression_results.Rdata",
    gene_data = "results/download_gene_data/gene_names.Rdata"
  output:
    factor_design_matrix = "results/limma_placentas/factor_design_matrix.csv",
    linear_models = "results/limma_placentas/fitted_model.Rdata",
    fdr_plot = "results/limma_placentas/fdr_plot.pdf",
    q_values_plot = "results/limma_placentas/q_value_distribution.pdf",
    q_values_plot_zoomed = "results/limma_placentas/q_value_distribution_zoomed.pdf",
    volcano_plots = "results/limma_placentas/volcano_plots.png",
    summary_csv = "results/limma_placentas/fold_change_summary.csv",
    summary_rds = "results/limma_placentas/fold_change_summary.rds",
    ranked_genes_upregulated = "results/limma_placentas/ranked_list_upregulated_genes.rds",
    ranked_genes_downregulated = "results/limma_placentas/ranked_list_downregulated_genes.rds"
  params:
    min_counts = MIN_COUNTS,
    fold_change_threshold = MIN_LOGFC,
    alpha = ALPHA
  script:
    "scripts/limma_placentas.R"

# Fit linear models for all non-placenta tissues
rule limma_non_placentas:
  input:
    expression_results = "results/tximeta/expression_results.Rdata",
    gene_data = "results/download_gene_data/gene_names.Rdata"
  output:
    factor_design_matrix = "results/limma_{tissue}/factor_design_matrix.csv",
    linear_models = "results/limma_{tissue}/fitted_model.Rdata",
    fdr_plot = "results/limma_{tissue}/fdr_plot.pdf",
    q_values_plot = "results/limma_{tissue}/q_value_distribution.pdf",
    q_values_plot_zoomed = "results/limma_{tissue}/q_value_distribution_zoomed.pdf",
    volcano_plots = "results/limma_{tissue}/volcano_plots.png",
    summary_csv = "results/limma_{tissue}/fold_change_summary.csv",
    summary_rds = "results/limma_{tissue}/fold_change_summary.rds",
    ranked_genes_upregulated = "results/limma_{tissue}/ranked_list_upregulated_genes.rds",
    ranked_genes_downregulated = "results/limma_{tissue}/ranked_list_downregulated_genes.rds"
  params:
    min_counts = MIN_COUNTS,
    fold_change_threshold = MIN_LOGFC,
    alpha = ALPHA,
    tissue = "{tissue}"
  wildcard_constraints:
    tissue = ".*_.*"
  script:
    "scripts/limma_non_placentas.R"
    
# Create variance stabilized counts for exploratory plotting
rule vsd:
  input:
    expression_results = 'results/tximeta/expression_results.Rdata'
  output:
    variance_stabilized_counts = 'results/tximeta/vsd.Rdata'
  params:
    min_counts = 9
  script:
    'scripts/vsd.R'

# Perform enrichment across all ranked sets of genes
rule enrichment_run:
  input:
    ranked_genes = "results/limma_{model}/ranked_list_{up_or_down}_genes.rds"
  output:
    raw_results = "results/gost_enrichment_run_{model}/gost_results_{up_or_down}.rds"
  params: 
    max_fdr = ALPHA
  script:
    "scripts/gost_enrichment_run.R"

# Format enrichment result to human readable tidy table
rule enrichment_format:
  input: "results/gost_enrichment_run_{model}/gost_results_{up_or_down}.rds"
  output: 
    go_long_format = "results/gost_enrichment_format/enrichment_{model}_{up_or_down}_long.csv",
    go_wide_format = "results/gost_enrichment_format/enrichment_{model}_{up_or_down}_wide.csv"
  script: 
    "scripts/gost_enrichment_format.R"
