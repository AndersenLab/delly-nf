
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

nextflow main.nf --debug

nextflow main.nf --sample_sheet=isotype_groups.tsv --species=c_elegans 

nextflow main.nf --sample_sheet=isotype_groups.tsv --bam_dir=/path/to/bams --reference=/path/to/reference.fa --ref_strain=reference_strain 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    (optional)
    --sample_sheet       TSV with isotype_ref_strain column (needs header)        (required)
    --minsize            The minimum size in bp to report for INDELs              (optional)
    --maxsize            The maximum size in bp to report for INDELs              (optional)
    
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     (required/optional)
    and / or
    --bam_dir            Path to folder containing bams                           (optional/required)
    --ref_strain         Name of strain to use a reference (matches genome ref)   (optional/required)
    --reference          Path to reference genome fasta                           (optional/required)
 

    HELP: https://andersenlab.org/dry-guide/latest/pipelines/pipeline-delly   
    ----------------------------------------------------------------------------------------------
"""

## Software Requirements

* The latest update requires Nextflow version 24+. On Rockfish, you can access this version by loading the `nf24_env` conda environment prior to running the pipeline command:

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
nextflow run -latest andersenlab/delly-nf --sample_sheet <path_to_sample_sheet> --bam_dir <path_to_bam_folder> --reference <path_to_reference> --ref_strain <reference_strain>
```


# Parameters

## -profile

There are three configuration profiles for this pipeline.

* `rockfish` - Used for running on Rockfish (default).
* `quest`    - Used for running on Quest.

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

A sample sheet produced by the concordance pipeline with a column specifying isotype reference strains with the column header "isotype_ref_strain".

## --species (required if --bam_dir, --reference, or --ref_strain not specified, otherwise optional)

Options: c_elegans, c_briggsae, or c_tropicalis

## --bam_dir (required if --species not specified, otherwise optional)

Path to the **folder** containing species strain bams.

## --reference (required if --species not specified, otherwise optional)

Path to the reference strain fasta file.

## --ref_strain (required if --species not specified, otherwise optional)

The name of the reference strain to call indels against

## --minsize (optional)

The minimum size cutoff for reporting an insertion or deletion (default: 50bp)

## --maxsize (optional)

The maximum size cutoff for reporting an insertion or deletion (default: 1000bp)

## -output-dir (optional)

__default__ - `results`

A directory in which to output results

# Output

```
└── results
    ├── workflow_software_versions.txt
    └── indels
        ├── AB1_indels.vcf.gz
        ├── AB1_indels.vcf.gz.tbi
        ├── AB4_indels.vcf.gz
        ├── AB4_indels.vcf.gz.tbi
        ├── MY23_indels.vcf.gz
        ├── MY23_indels.vcf.gz.tbi
        ...
```

# Relevant Docker Images

* `dellytools/delly` ([link](https://hub.docker.com/r/dellytools/delly))
* `quay.io/biocontainers/bcftools` ([link](https://quay.io/biocontainers/bcftools))


*Note: If running on Rockfish, make sure to properly set up Nextflow prior to running workflow ([Nextflow configuration](http://andersenlab.org/dry-guide/latest/rockfish/rf-nextflow/#configuring_nextflow)).*
