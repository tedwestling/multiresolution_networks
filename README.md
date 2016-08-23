# multiresolution_networks
Replication files for multiresolution networks paper

These files contain all code necessary to replicate the results presented in “Multiresolution network models” (Fosdick, McCormick, Murphy, Ng, and Westling). 

Prior to running the code in scripts/indian_village_estimation.R, the user must first:

1) download the Karnataka village dataset (as a zip file) from https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21538

2) unzip the file

3) place the directories ‘Data/1. Network Data’ and ‘Data/2. Demographics and Outcomes’ from the unzipped file in data/indian_village_raw

Once the data is in place, the driver script for running all the code can be found in code/scripts/indian_village_estimation.R. This code will choose the number of blocks, run the MCMC sampling algorithm, post-process the samples, and make all plots that appear in the paper. Both the choice of blocks and the MCMC sampling take a long time to run. The data ultimately created by the code is also already present in the data/results directory, so it is not necessary to run the entire block choice or MCMC.

Note that you need the April 2016 version of the smacof package, otherwise you'll get an error when running the postprocessing script.