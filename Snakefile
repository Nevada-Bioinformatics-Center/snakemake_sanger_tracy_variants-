import glob

configfile: "config.yaml"

inputdirectory=config["input_data"]
SAMPLES, =glob_wildcards(inputdirectory+"/{sample}.ab1", followlinks=True)

wildcard_constraints:
    sample = "[A-Za-z0-9\-]+"

##### target rules #####
rule all:
    input: 
       expand("variants/{sample}.fasta", sample=SAMPLES),
       expand("variants/{sample}.bcf", sample=SAMPLES),
       #expand("variants/{sample}.vcf.gz.tbi", sample=SAMPLES),
       expand("variants/{sample}.snpEff.vcf", sample=SAMPLES),
       expand("variants/{sample}.snpEff.html", sample=SAMPLES),
       expand("variants/{sample}.snpEff.vcf.gz.tbi", sample=SAMPLES),
       expand("variants/{sample}.highQ.snpEff.vcf", sample=SAMPLES),
       expand("variants/{sample}.highQ.snpEff.html", sample=SAMPLES),
       expand("variants/{sample}.highQ.snpEff.vcf.gz.tbi", sample=SAMPLES),
       "variants/all.vcf",
       "variants/all.snpEff.vcf.gz",
       "variants/all.snpEff.vcf.gz.tbi",
       "variants/all.highQ.snpEff.vcf.gz",
       "variants/all.highQ.snpEff.vcf.gz.tbi",



rule tracy_call_varaints:
###ref must be indexed with tracy prior to running pipeline
    input:
       ref=config["ref"],
       #ref="ref/chr1_chr5.fasta.gz",
       samplefile=f"{config['input_data']}/{{sample}}.ab1"
    output:
        "variants/{sample}.bcf"
    params:
        "variants/{sample}"
    log:
        "logs/tracy/variants/{sample}.log",
    conda:
        "envs/tracy.yaml"
    shell:
        "tracy decompose -v -o {params} -r {input.ref} {input.samplefile} > {log} 2>&1"

rule tracy_fasta:
    input:
       samplefile=f"{config['input_data']}/{{sample}}.ab1"
    output:
        "variants/{sample}.fasta"
    log:
        "logs/tracy/fasta/{sample}.log"
    conda:
        "envs/tracy.yaml"
    shell:
        "tracy basecall -f fasta -o {output} {input.samplefile} > {log} 2>&1"

rule bcf_to_vcf:
    input:
        bcf="variants/{sample}.bcf"
    output:
        vcf="variants/{sample}.temp.vcf"
    params:
        ""  # optional parameters for bcftools view (except -o)
    log:
        "logs/bcftools/{sample}.log"
    wrapper:
        "v1.3.2/bio/bcftools/view"

rule vcf_add_samplename:
    input:
        "variants/{sample}.temp.vcf"
    output:
        "variants/{sample}.vcf"
    params:
        "{sample}"
    log:
        "logs/sed/{sample}.log"
    conda:
        "envs/sed.yaml"
    shell:
        """
        sed 's/sample/{params}/g' {input} > {output} 2> {log}
        rm {input}
        """

rule bgzip_files:
    input:
        "variants/{sample}.vcf"
    output:
        "variants/{sample}.vcf.gz"
    log:
        "logs/bgzip/{sample}.log"
    params:
        extra="", # optional
    threads: 1
    conda:
        "envs/bgzip.yaml"
    shell:
        "bgzip {input} 2> {log}"

rule index_files:
    input:
        "variants/{sample}.vcf.gz"
    output:
        "variants/{sample}.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{sample}.log"
    wrapper:
        "v1.3.2/bio/tabix"

rule snpeff_annotate_persample:
    input:
        "variants/{sample}.vcf.gz"
    output:
        vcf="variants/{sample}.snpEff.vcf",
        html="variants/{sample}.snpEff.html",
        txt="variants/{sample}.snpEff.txt",
    params:
        dbname=config["snpeffname"]
    log:
        "logs/snpeff/{sample}.log"
    conda:
        "envs/snpeff.yaml"
    shell:
        #snpEff Oryza_sativa {input} | bcftools view -Oz > {output} 2> {log}
        """
        snpEff {params.dbname} {input} > {output.vcf} 2> {log} && mv snpEff_genes.txt {output.txt} && mv snpEff_summary.html {output.html}
        """

rule snpeff_annotate_all:
    input:
        "variants/all.vcf"
    output:
        vcf="variants/all.snpEff.vcf",
        html="variants/all.snpEff.html",
        txt="variants/all.snpEff.txt",
    params:
        #cf=config["snpeffconfig"],
        dbname=config["snpeffname"]
    log:
        "logs/snpeff/all.log"
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        snpEff {params.dbname} {input} > {output.vcf} 2> {log} && mv snpEff_genes.txt {output.txt} && mv snpEff_summary.html {output.html}
        """

rule highq_vcf_persample:
    input:
        "variants/{sample}.vcf.gz"
    output:
        "variants/{sample}.highQ.vcf.gz"
    params:
        extra="-f PASS"
    log:
        "logs/bcftools/filter_{sample}.log"
    wrapper:
        "v1.3.2/bio/bcftools/view"

rule index_files_highq:
    input:
        "variants/{sample}.highQ.vcf.gz"
    output:
        "variants/{sample}.highQ.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{sample}.log"
    wrapper:
        "v1.3.2/bio/tabix"

rule snpeff_annotate_persample_highq:
    input:
        "variants/{sample}.highQ.vcf.gz"
    output:
        vcf="variants/{sample}.highQ.snpEff.vcf",
        html="variants/{sample}.highQ.snpEff.html",
        txt="variants/{sample}.highQ.snpEff.txt",
    params:
        dbname=config["snpeffname"]
    log:
        "logs/snpeff/{sample}_highq.log"
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        snpEff {params.dbname} {input} > {output.vcf} 2> {log} && mv snpEff_genes.txt {output.txt} && mv snpEff_summary.html {output.html}
        """

rule merge_vcf:
    input:
        calls=expand("variants/{sample}.vcf.gz", sample=SAMPLES),
        tbi=expand("variants/{sample}.vcf.gz.tbi", sample=SAMPLES)
    output:
        vcf="variants/all.vcf"
    log:
        "logs/bcftools/merge_all.log"
    wrapper:
        "v1.3.2/bio/bcftools/merge"

rule merge_vcf_highQ:
    input:
        calls=expand("variants/{sample}.highQ.vcf.gz", sample=SAMPLES),
        tbi=expand("variants/{sample}.highQ.vcf.gz.tbi", sample=SAMPLES)
    output:
        vcf="variants/all.highQ.vcf"
    log:
        "logs/bcftools/merge_all_highQ.log"
    wrapper:
        "v1.3.2/bio/bcftools/merge"

rule snpeff_annotate_all_highq:
    input:
        "variants/all.highQ.vcf"
    output:
        vcf="variants/all.highQ.snpEff.vcf",
        html="variants/all.highQ.snpEff.html",
        txt="variants/all.highQ.snpEff.txt",
    params:
        dbname=config["snpeffname"],
    log:
        "logs/snpeff/all_highq.log",
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        snpEff {params.dbname} {input} > {output.vcf} 2> {log} && mv snpEff_genes.txt {output.txt} && mv snpEff_summary.html {output.html}
        """


rule bgzip_files_all:
    input:
        "variants/all.snpEff.vcf"
    output:
        "variants/all.snpEff.vcf.gz"
    log:
        "logs/bgzip/all.log"
    params:
        extra="", # optional
    threads: 1
    conda:
        "envs/bgzip.yaml"
    shell:
        "bgzip {input} 2> {log}"

rule index_files_all:
    input:
        "variants/all.snpEff.vcf.gz"
    output:
        "variants/all.snpEff.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/all.log"
    wrapper:
        "v1.3.2/bio/tabix"

#rule bgzip_files_all_highq:
#    input:
#        "variants/all.snpEff.highQ.vcf"
#    output:
#        "variants/all.snpEff.highQ.vcf.gz"
#    log:
#        "logs/bgzip/all_snpeff_highq.log"
#    params:
#        extra="", # optional
#    threads: 1
#    conda:
#        "envs/bgzip.yaml"
#    shell:
#        "bgzip {input} 2> {log}"

#rule index_files_all_highq:
#    input:
#        "variants/all.snpEff.highQ.vcf.gz"
#    output:
#        "variants/all.snpEff.highQ.vcf.gz.tbi"
#    params:
#        "-p vcf"
#    log:
#        "logs/tabix/all_snpeff_highq.log"
#    wrapper:
#        "v1.3.2/bio/tabix"

rule bgzip_files_all_highq:
    input:
        "variants/all.highQ.snpEff.vcf"
    output:
        "variants/all.highQ.snpEff.vcf.gz"
    log:
        "logs/bgzip/all_snpeff_highq.log"
    params:
        extra="", # optional
    threads: 1
    conda:
        "envs/bgzip.yaml"
    shell:
        "bgzip {input} 2> {log}"


rule index_files_all_highq:
    input:
        "variants/all.highQ.snpEff.vcf.gz"
    output:
        "variants/all.highQ.snpEff.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/all_snpeff_highq.log"
    wrapper:
        "v1.3.2/bio/tabix"
