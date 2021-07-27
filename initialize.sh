#!/bin/bash
Projecttitle="UC-RNAseq"
User="Eike"

echo $'# README\n\n## Table of Contents\n\n| File Name | Description |\n| -- | -- |' > Notebook_${User}/README.md
rm -r 01_berrycounts 02_berryreplace 03_berryhash Notebook_Eike Notebook_Florian
rm 00_RawData/results/berries.txt 
mkdir Notebook_${User}
sed -i -n '/#/p' 00_Files.md
sed -i -n '/#/p' 01_Background.md
echo '# ${!Projecttitle}' > README.md