import pandas as pd

configfile : "config.yaml"
 
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

conditions = pd.read_table("samples.tsv").set_index("condition", drop=False)

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml") 
 


 
include: "rules/00_download_data.skm"

#################################### Mapping and Quantification ################################
#
# In this module, we are declaring four rules that are designed to map all the reads to the  
# genome (hisat2) and count the reads that map to each gene (featureCounts). 
#
#########################################################################################    
 

    
rule hisat2_Genome_index:  #This is a rule and represent the first step of mapping the reads with hisat (indexing the genome)
    input:
        "Genome/" + config["assembly"] + ".fa"
    output:
        "Genome/Index/" + config["assembly"] + ".1.ht2"
    threads: 7
    conda:
        "envs/core.yaml"
    log:
        "logs/hisat2_Genome_index.log"
    shell:
        "hisat2-build -p {threads} {input} Genome/Index/" + config["assembly"]  + " 2> {log}"




def sample_to_unit(wildcards):
    return units.loc[(wildcards.sample, "1" ) , ("fq1", "fq2") ].dropna() # We are not yet supporting for lanes	

#def get_fastq(wildcards):
#    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


	
if str2bool(config["paired_end"])==False:
        
    rule hisat2_to_Genome:
        input:
            fastq = sample_to_unit,
            genome = "Genome/Index/" + config["assembly"] + ".1.ht2"
        output:
            temp("hisat2/{sample}.sam")
        threads: 6
        log:
            "logs/hisat2_{sample}.log"       
        conda:
            "envs/core.yaml"
        shell:
            "hisat2 -p {threads} -U {input.fastq} -x  Genome/Index/" + config["assembly"] +  "  > {output}  2> {log} "
            
elif str2bool(config["paired_end"])==True:
    
    rule hisat2_to_Genome:
        input:
            fastq = sample_to_unit,
            genome = "Genome/Index/" + config["assembly"] + ".1.ht2"
        output:
            temp("hisat2/{sample}.sam")
        threads: 6
        log:
            "logs/hisat2_{sample}.log"    
        conda:
            "envs/core.yaml"
        shell:
            "hisat2 -p {threads} -1 {input.fastq[0]} -2 {input.fastq[1]} -x  Genome/Index/" + config["assembly"] +  "  > {output}  2> {log} "


rule samTobam:
    input:
        "hisat2/{sample}.sam"
    output:
        "hisat2/{sample}.sorted.bam"
    conda:
        "envs/core.yaml"
    shell:
        "samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output} "
        
rule bamstats:
    input:
        "hisat2/{sample}.sorted.bam"
    output:
        stats_txt = "QC/{sample}/{sample}.stats",
        stats_html = "QC/{sample}/{sample}.plots.html"
    params:
        "QC/{sample}/{sample}.plots"
    conda:
        "envs/core.yaml"
    shell:
        "samtools stats {input} > {output.stats_txt} && plot-bamstats -p {params} {output.stats_txt}"

      
########

#rule featureCounts:
#    input:
        #gtf = "gffcompare/extended_ref_annotation.gtf",
#        gtf = "Gene_annotation/" + config["assembly"] + ".ensGene.gtf",
#        bam = expand("hisat2/{sample}.sorted.bam", sample=SAMPLES)
#    output:
#        "featureCounts/total_samples.gene_count.txt"
#    threads: 1
#    conda:
#        "envs/core.yaml"
#    log:
#        "logs/featureCounts.total.log"
#    shell:
#        "featureCounts -a {input.gtf} -o {output} {input.bam} 2> {log}"

 
 
rule featureCounts:
    input:
        gtf = "Gene_annotation/" + config["assembly"] + ".ensGene.gtf",
        bam = "hisat2/{sample}.sorted.bam"
    output:
        "featureCounts/{sample}.gene_count.txt"
    threads: 1
    conda:
        "envs/core.yaml"
    log:
        "logs/featureCounts.{sample}.log"
    shell:
        "featureCounts -a {input.gtf} -o {output} {input.bam} 2> {log}"
 
############# Downstream analysis #############
#
# Everything below corresponds to workflows to perform different anlyses to get meaningful 
# quantitative data. On rules/ folder you can see the different snakemake modules (.skm files)
# which are `included` to be connected with the previous rules that are explicit on this
# current script. The `include` statement allows the integration of the .skm files. Notice 
# that all these snakemake scripts work under python, thus any python syntax can be used.
# 
###############################################    

#####  DGA

include: "rules/diffexp.smk"

rule run_DGA:
   input:
       expand(["results/diffexp/{contrast}.diffexp.tsv",
               "results/diffexp/{contrast}.ma-plot.svg"],
              contrast=config["diffexp"]["contrasts"])

######






include: "rules/Pseudoalignment.skm"    
     
rule run_salmon:
    input:
        expand( 'salmon/{sample}/quant.sf', sample=SAMPLES)
    
rule genecount:
    input:
        "featureCounts/total_samples.gene_count.txt", 
        expand( 'salmon/{sample}/quant.sf', sample=SAMPLES)     
    
#include: "rules/01_stringtie.skm"    
#include: "rules/02_bridge.skm"  
#include: "rules/03_whippet_quant.skm"
include: "rules/03.1_whippet_quant.skm"

#rule get_whippet_quant:    #This is a calling point to run all whippet analysis
#    input:
#        expand("Whippet/Quant/{sample}.psi.gz", sample=SAMPLES)
    
#include: "rules/04_whippet_delta.skm"
include: "rules/04.1_whippet_delta.skm" 

    

