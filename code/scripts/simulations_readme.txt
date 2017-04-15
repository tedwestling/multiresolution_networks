LS-SBM SIMULATIONS README

Simulations were run on a cluster with slurm task management. The simulation directory should contain:
- the file run_one_sim.R
- the file header.R
- a subdirectory code/functions/ containing all the project function files
- empty subdirectories outs/ and results/
- and a file run_sims.sbatch with the following code:

########################################################
#!/bin/bash
#SBATCH --job-name multires_sim 
#SBATCH --partition short
#SBATCH --ntasks 1  
#SBATCH --cpus-per-task=1
#SBATCH --time 0-12:00      
#SBATCH --mem-per-cpu=2000     # Memory limit for each tasks (in MB)
#SBATCH -o outs/multires_sim_%a.out    # File to which STDOUT will be written
#SBATCH -e outs/multires_sim_%a.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tgwest@uw.edu # Email to which notifications will be sent

R --vanilla --no-restore --no-save --args n ${1} < run_one_sim.R > outs/multires_sim_"${SLURM_ARRAY_TASK_ID}".Rout
########################################################

The simulations are then run by typing the following in to the command line:

sbatch --array=1-1000 run_sims.sbatch 300