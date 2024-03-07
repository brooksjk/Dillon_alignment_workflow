# Ancient/Modern Alignment Workflow
## Files included 
- *slurm/config.yaml*: config file for HPC architecture and slurm capability
- *Snakefile*: the pipeline
- *config.yaml*: environmental variables for pipeline
- *initiator.sh*: setps up directory and launches snakemake_submitter.sh
- *snakemake_submitter.sh*: initiates conda environment and submits snakemake job
  
## Getting Started
To ensure that the workflow functions correctly, you will need to make some minor adjustments. Specifically, you will need to modify the config.yaml file and the snakemake_submitter files. The changes you need to make are indicated within angle brackets (< >) in these files.

If you would like to ensure that the workflow is running correctly after you have made your modifications a dry run can be conducted using the following commands on Secretariat:
```
#load in conda environments and activate snakemake
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
conda activate snakemake

#Snakemake dry run command
#Needs to be done in directory containing Snakefile
snakemake -n 
```
If done correctly you should see the following:
- Workflow Summary: Snakemake should display a summary of the workflow, including the total number of jobs that are ready to be executed.
- Job Details: For each job, Snakemake will list the rule name, input files, output files, and the shell command that will be executed. Ensure that these details align with your expectations.
- No Errors: There should be no error messages or warnings indicating missing input files, syntax errors, or other issues.

If all these conditions are met, your workflow is likely set up correctly and ready for execution.

## Running the workflow
In order to submit the workflow to run you will need to move into your working directory and submit the following command on Secretariat: `sbatch initiator.sh`

If some jobs finish, but there are a few incomplete jobs you can modify the Snakemake statement in the snakemake_submitter.sh file to include the --rerun-incomplete option as such:
```
snakemake \
--rerun-incomplete \
-s Snakefile \
--profile slurm \
--latency-wait 150 \
-p \
```



