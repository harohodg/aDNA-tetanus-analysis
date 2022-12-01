This folder includes the main scripts which were run on Compute Canada (now the [Digital Research Alliance of Canada](https://alliancecan.ca/en)) high performance computing resources. 

`WARNING` : Most of these scripts assume they were launched on a Compute Canada system with `sbatch` and won't work properly if run directly or on another system.


Most programs have two scripts. One name `program_job.sh` and one named `launch_program_jobs.sh`. The former runs a single instance of the program with a specific set of parameters, the latter takes an `input folder`, `output folder`, and `set of parameters` and launches a collection of jobs with the same set of parameters for each input file found in the input folder. 

All scripts should be run via `sbatch` or with an interactive job which can be launched with `launch_interactive_job.sh`
Ideally all of these scripts would be bundled together with apptainer containers and run with something like Nextflow


Most of these scripts should run on any of the major systems including [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en), [Cedar](https://docs.alliancecan.ca/wiki/Cedar), [Graham](https://docs.alliancecan.ca/wiki/Graham), and [Narval](https://docs.alliancecan.ca/wiki/Narval/en) with the following exceptions.
- any script for downloading SRA datasets should only be run on Cedar compute nodes to avoid overwhelming the login nodes.
- The MEGAHIT scripts were coded to work on and were tested on [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en). They try to determine how much RAM will be needed and if it's above 80G they request an entire node. 


The apptainer `*.sif` files for mapDamage and pyDamage were built using `apptainer build mapDamage.sif docker://quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0` and `apptainer build pyDamage.sif docker://quay.io/biocontainers/pydamage:0.70--pyhdfd78af_1`



