date
echo "Setting up R 4.4.0 & Python 3.11.5 Environments"
# load environment and programming languages
module load StdEnv/2023
module load gcc/12.3
module load gdal/3.7.2
module load r/4.4.0
module load python/3.11.5

pwd 

# If getting errors when booting up R (e.g. LC_XXX failed, using "C"), 
# run the following or add to .bash_profile
# export LANG="en_US.UTF-8" 