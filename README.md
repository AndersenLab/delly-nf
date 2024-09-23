![Build Docker (env/Dockerfile)](https://github.com/AndersenLab/delly-nf/workflows/Build%20Docker%20(env/Dockerfile)/badge.svg)


# delly-nf

This pipeline performs INDEL calling on isotype strains versus the reference strain. The VCFs output from this pipeline are used within the lab and also released to the world via CaeNDR.

# Pipeline overview

"""

      ______   _______  _        _                   _        _______ 
     (  __  \ (  ____ \( \      ( \   |\     /|     ( (    /|(  ____ \
     | (  \  )| (    \/| (      | (   ( \   / )     |  \  ( || (    \/
     | |   ) || (__    | |      | |    \ (_) /_____ |   \ | || (__    
     | |   | ||  __)   | |      | |     \   /(_____)| (\ \) ||  __)   
     | |   ) || (      | |      | |      ) (        | | \   || (      
     | (__/  )| (____/\| (____/\| (____/\| |        | )  \  || )      
     (______/ (_______/(_______/(_______/\_/        |/    )_)|/       
                                                                                                                                                                                                                
nextflow main.nf --help

nextflow main.nf -profile rockfish --debug

nextflow main.nf -profile rockfish --sample_sheet=isotype_groups.tsv --species=c_elegans 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    (optional)
    --sample_sheet       TSV with column isotype (needs header)                   (required)
    --masking            BED file containing regions to skip during indel calling (optional)
    --minsize            The minimum size in bp to report for INDELs              (optional, default: 50bp)
    --maxsize            The maximum size in bp to report for INDELs              (optional, default: 1000bp)
    --output             Output folder name (optional)                            (optional)
    
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     (required/optional)
    or
    --bam_dir            Path to folder containing bams                           (optional/required)
    --reference          Path to reference genome fasta                           (optional/required)
 

    HELP: https://andersenlab.org/dry-guide/latest/pipelines/pipeline-delly   
    ----------------------------------------------------------------------------------------------
"""

## Software Requirements

* The latest update requires Nextflow version 23+. On Rockfish, you can access this version by loading the `nf23_env` conda environment prior to running the pipeline command:

```
module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env
```

# Usage

*Note: if you are having issues running Nextflow or need reminders, check out the [Nextflow](http://andersenlab.org/dry-guide/latest/rockfish/rf-nextflow/) page.*

## Testing on Rockfish

*This command uses a test dataset*

```
nextflow run -latest andersenlab/delly-nf --debug
```

## Running on Rockfish

You should run this in a screen or tmux session.

```
nextflow run -latest andersenlab/delly-nf --sample_sheet <path_to_sample_sheet> --species <species>
```

or

```
nextflow run -latest andersenlab/delly-nf --sample_sheet <path_to_sample_sheet> --bam_dir <path_to_bam_folder> --reference <path_to_reference>
```


# Parameters

## -profile

There are three configuration profiles for this pipeline.

* `rockfish` - Used for running on Rockfish (default).
* `quest`    - Used for running on Quest.
* `local`    - Used for local development.

>[!Note]
>If you forget to add a `-profile`, the `rockfish` profile will be chosen as default

## --debug

You should use `--debug` for testing/debugging purposes. This will run the debug test set (located in the `test_data` folder).

For example:

```
nextflow run -latest andersenlab/delly-nf --debug
```

Using `--debug` will automatically set the sample sheet to `test_data/sample_sheet.tsv`

## --sample_sheet

A sample sheet produced by the concordance pipeline with a column specifying isotype reference strains with the column header "isotype".

## --species (required if --bam_dir and --reference not specified, otherwise optional)

Options: c_elegans, c_briggsae, or c_tropicalis

## --bam_dir (required if --species not specified, otherwise optional)

Path to the **folder** containing species strain bams.

## --reference (required if --species not specified, otherwise optional)

Path to the reference strain fasta file.

## --masking (optional)

Path to bed file containing regions to be skipped during INDEL calling. For C. elegans this defaults to HVR calls in test_data/c_elegans_mask.bed.

## --minsize (optional)

The minimum size cutoff for reporting an insertion or deletion (default: 50bp)

## --maxsize (optional)

The maximum size cutoff for reporting an insertion or deletion (default: 1000bp)

## --output (optional)

__default__ - `delly-YYYYMMDD`

A directory in which to output results. If you have set `--debug`, the default output directory will be `delly-YYYYMMDD-debug`.

# Output

```
└── ANNOTATE_VCF
    ├── AB1_indels_filtered.vcf.gz
    ├── AB1_indels_filtered.vcf.gz.tbi
    ├── AB1_indels_unfiltered.vcf.gz
    ├── AB1_indels_unfiltered.vcf.gz.tbi
    ├── MY23_indels_filtered.vcf.gz
    └── MY23_indels_filtered.vcf.gz.tbi
    ├── MY23_indels_unfiltered.vcf.gz
    └── MY23_indels_unfiltered.vcf.gz.tbi

```

# Relevant Docker Images

* `andersenlab/delly` ([link](https://hub.docker.com/r/andersenlab/delly)): Docker image is created within this pipeline using GitHub actions. Whenever a change is made to `env/Dockerfile` or `.github/workflows/build_docker.yml` GitHub actions will create a new docker image and push if successful
* `andersenlab/annotation` ([link](https://hub.docker.com/r/andersenlab/annotation)): Docker image is created within this pipeline using GitHub actions. Whenever a change is made to `env/annotation.Dockerfile` or `.github/workflows/build_docker.yml` GitHub actions will create a new docker image and push if successful


Make sure that you add the following code to your `~/.bash_profile`. This line makes sure that any singularity images you download will go to a shared location on `/vast/eande106` for other users to take advantage of (without them also having to download the same image).

```
# add singularity cache
export SINGULARITY_CACHEDIR='/vast/eande106/singularity/'
```
