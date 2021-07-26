#!/usr/bin/env bash
#purpose of the script
#This is a shell script to replicate the whole step if need be. 
#in this script, necessary environments should be set. A complete conda environment.yml should exist in the main project directory.
#this script saves a result to the results/ directory. Data ok to be on the repo is exported to results/. Anything to large or private should go to local/. Anything safe to delete after running is going to temp/. Anything which is a major outcome of the project is noted down in 00_Files.md
conda activate berry_env
mkdir results/ local/ temp/ #whatever you need.
#your commands start here
#Rscript script.R
#...
