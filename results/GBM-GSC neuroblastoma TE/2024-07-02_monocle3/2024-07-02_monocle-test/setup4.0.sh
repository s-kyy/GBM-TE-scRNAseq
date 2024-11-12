date
echo "Setting up R 4.0.2 & Python 3.8.2 Environments"
# load environment and programming languages
module load StdEnv/2020; module load nixpkgs/16.09 gcc/7.3.0
module load r/4.0.2  
module load python/3.8.2 

pwd 

# If getting errors when booting up R (e.g. LC_XXX failed, using "C"), 
# run the following or add to .bash_profile
# export LANG="en_US.UTF-8" 
