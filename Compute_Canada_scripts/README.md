This folder includes the main analysis scripts which were run on Compute Canada (now the [Digital Research Alliance of Canada](https://alliancecan.ca/en)) high performance computing resources. 

Most programs have two scripts. One name `program_job.sh` and one named `launch_program_jobs.sh`. The former runs a single instance of the program with a specific set of parameters, the latter takes an `input folder`, `output folder`, and `set of parameters` and launches a collection of jobs with the same set of parameters for each input file found in the input folder. 

`WARNING` : Most of these scripts assume they were launched on a Compute Canada system with `sbatch` and won't work properly if run directly or on another system. 

Most of these scripts should run on any of the major systems including [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en), [Cedar](https://docs.alliancecan.ca/wiki/Cedar), [Graham](https://docs.alliancecan.ca/wiki/Graham), and [Narval](https://docs.alliancecan.ca/wiki/Narval/en) with the following exceptions.
- `download_SRA_datasets.sh`, `fasterq-dump_job.sh`, `prefetch_SRA_dataset.sh`, and `prefetch_SRA_datasets.sh` should only be run on Cedar compute nodes to avoid overwhelming the login nodes.
- `launch_megahit_jobs.sh` and `launch_megahit_coassembly_jobs.sh` were coded to work on and were tested on [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en). They try to determine how much RAM will be needed and if it's above 80G they request an entire node. 

