#!/bin/bash
mkdir Notebook_User
echo $'# README\n\n## Table of Contents\n\n| File Name | Description |\n| -- | -- |' > Notebook_User/README.md
rm -r 01_berrycounts 02_berryreplace 03_berryhash Notebook_Eike Notebook_Florian
rm 00_RawData/results/berries.txt 
sed -i -n '/#/p' 00_Files.md
sed -i -n '/#/p' 01_Background.md
echo '# Your new project' > README.md