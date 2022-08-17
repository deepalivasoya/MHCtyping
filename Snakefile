configfile: "config.yaml"
# print(config['samples'][0])

rule main:
  input: expand("results/{sample}/{sample}.{pcr}.clusters.blast.log", sample=config["samples"], pcr=config["primers"].keys())

#Rule for quality trimming sequencing reads
rule trim:
 output:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
  singles = "results/{sample}/{sample}.trimmed_singles.fastq",
  log = "results/{sample}/{sample}.trimming_by_sickle.log"
 input:
  read1 =lambda wildcards: "fastq/{}_1.fastq.gz".format(config['reads'][wildcards.sample]),
  read2 =lambda wildcards: "fastq/{}_2.fastq.gz".format(config['reads'][wildcards.sample]),
 params:
  qual_threshold = config["read_quality_threshold"],
  len_threshold = config["read_trim_length_threshold"],
  folder = "results/{sample}"
 shell:
  r"""sickle pe -f {input.read1} -r {input.read2} -o {output.read1} -p {output.read2} -s {output.singles} -t sanger -q {params.qual_threshold} -l {params.len_threshold} > {output.log}
    fastqc -f fastq -o {params.folder} {input.read1}
    fastqc -f fastq -o {params.folder} {input.read2}
  """

#Rule to extend overlapping paired end reads
rule overlap:
 output: "results/{sample}/{sample}.overlap_by_flash.log"
 input:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
 params:
  min_overlap = config["minimum_overlap"],
  max_overlap = config["maximum_overlap"],
  directory = "results/{sample}",
  prefix = "{sample}"
 shell:
  r"""flash -m {params.min_overlap} -M {params.max_overlap} -O -o {params.prefix} -d {params.directory} --threads=1 {input.read1} {input.read2} > {output}
      gzip {input.read1}
      gzip {input.read2}
      gzip results/{params.prefix}/{params.prefix}.notCombined_1.fastq
      gzip results/{params.prefix}/{params.prefix}.notCombined_2.fastq
      gzip results/{params.prefix}/{params.prefix}.trimmed_singles.fastq
  """

#Rule to segragate extended reads using primer sequencing and the cluster 100% identical sequences
rule sort:
 output: "results/{sample}/{sample}.{pcr}.sort.stats.csv"
 input: "results/{sample}/{sample}.overlap_by_flash.log"
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  pcr = "{pcr}",
  primers = lambda wildcards: config['primers'][wildcards.pcr]
 shell:
  "perl scripts/sortPrimers.pl --sample={params.prefix} --work_dir={params.directory} --primers={params.pcr} --primer_seq={params.primers}"

#Rule to check each cluster (top 1000 clusters only) for chimera, 1bp variants, >9bp length difference, ambiguous calls etc. 
rule filtering:
 output: "results/{sample}/{sample}.{pcr}.clusters.stats.csv"
 input: "results/{sample}/{sample}.{pcr}.sort.stats.csv"
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  pcr = "{pcr}",
  length = lambda wildcards: config['amplicon_size'][wildcards.pcr],
  no_of_clusters = config["no_of_clusters"]
 shell:
  "perl scripts/filterSequences.pl --sample={params.prefix} --work_dir={params.directory} --primers={params.pcr} --ampliconSize={params.length}"

#Rule to run blast 
rule blast:
 output: "results/{sample}/{sample}.{pcr}.clusters.blast"
 input: 
  log = "results/{sample}/{sample}.{pcr}.clusters.stats.csv"
 params:
  reference = lambda wildcards: config['database'][wildcards.pcr],
  query = "results/{sample}/{sample}.{pcr}.clusters.fasta"
 shell:
  "blastn -db {params.reference} -query {params.query} -outfmt '6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore' -out {output}"
  
#Rule for analying the blast result and categoried the result
rule analyse_blast:
 output: "results/{sample}/{sample}.{pcr}.clusters.blast.stats.csv"
 input: "results/{sample}/{sample}.{pcr}.clusters.blast"
 params:
  prefix = "{sample}",
  pcr = "{pcr}",
  directory = "results/{sample}",
  filtered = "results/{sample}/{sample}.{pcr}.clusters.details.txt",
  clusters_fasta = "results/{sample}/{sample}.{pcr}.clusters.fasta",
  cutoff = lambda wildcards: config['read_count_percent_threshold'][wildcards.pcr],
  fc = lambda wildcards: config['fold_change_discarded'][wildcards.pcr]
 shell:
  "perl scripts/reportMapping.pl --sample={params.prefix} --work_dir={params.directory} --primers={params.pcr} --blast_file={input} --cluster_details={params.filtered} --cluster_fasta={params.clusters_fasta} --cutoff={params.cutoff} --fc={params.fc}"
  



