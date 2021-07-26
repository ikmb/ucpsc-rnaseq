## Replace the word "berry"

*20.07.2021*

user@com123:~/Desktop/uc-rnaseq/02_berryreplace  

I need to replace the occurence of the word berry in the names of the berries to feed it to a program allergic to berries.

    bash runReplaceBerries.sh  
	cat local/replaced\_berries_result.txt  
	berry,number
	rasp,4
	straw,14
	blue,6

Result:  
Except for the header line, all instances of "berry" were successfully removed. I saved the result to local/ to demonstrate how we would deal with files which can not be part of the repo due to size or privacy. It's important that the path where the script was run is given in the wrapper shell script of the step. This way, hidden dependency problems occuring later can be investigated more easily. In the notebook, the path has to be given anyway. The new file's path is noted down in the 00_Files.md.
