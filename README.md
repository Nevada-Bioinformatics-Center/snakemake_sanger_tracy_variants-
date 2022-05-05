# Snakemake workflow: snakemake_sanger_tracy_variants

This workflow performs SNP calling of .ab1 files from sanger sequencing to a reference genome. The resulting VCF files are then ran through snpEff to predict
variant effects and create summary .html and .txt files. The data is also filtered for high quality mutants and result files are labeled as such.

## Authors

* Hans Vasquez-Gross (@hansvg)

## Usage

### Simple

#### Step 1: Install workflow

clone this workflow to your local computer

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

The `ref` directory variable must contain preindexed reference genome from Tracy

    tracy index -o hg38.fa.fm9 hg38.fa.gz
    samtools faidx hg38.fa.gz

The `input_data` variable must be a directory with all ab1 files for a project.

The `snpeffname` must be the snpEff database name pre-downloaded before running the pipeline.

    snpEff -download Oryza_sativa

You can create a conda environment with all the necessary tools as follows:

    mamba create -c conda-forge -c bioconda -n sanger_tools tracy snpeff samtools bcftools tabix


#### Step 3: Execute workflow

Test your configuration by performing a dry-run via 

    snakemake --use-conda --jobs 1 -prn

if no errors are presented, run the pipeline with the following command

    snakemake --use-conda --jobs 1 -pr

