# Docu-Template

This is a demo repo as an example for how scientific documentation in data science projects can function. To use this as a template for your data science project, select "Use this template" and open your repo. To remove all demo content from the repo, run `./initialize.sh` in your new repo's main folder.

## Layout description

The repo has three main areas:

* Main directory
* Researchers' notebooks
* Subdirectories of each independent analytis step

## Main directory

Latest when external collaborators join the project or it is to be published, the landing page of the repo should look welcoming and helpful. A person not knowing the project should be able to understand and reproduce the scientific work presented. The two main design elements for this are obvious filename and directory structure, and an introductory Readme.md.  

### Readme.md

The Readme.md includes a description of the project and points at the points where to get started, if this is the first contact with the repo:

* 00\_Files.md, 01_Background.md
* Notebooks
* Analysis directories
* Analysis flowchart?
* For publication, an as-short-as-possible but complete step by step guide to reproduce the results.

### 00_Files.md

This describes all files of importance for this project, namely:

* Input files, structured by type of data, source of data, experimental method as one sees fit. 
* Processed data which is not saved in the repo due to privacy or size, linking to databases, temp-folders, who to ask for access.
* Final Files of Importance. Everything that is considered an informative result and definitely all items that are included in the paper or supplement.

### 00a_Metadata.md

All information about files that can not easily be formatted into a table.

* Method details before writing the method section
* Method details so specific that they need not be published

A table of individual mice and their weight used as covariate for a sequencing experiment should be saved as a file.

### 01_Background.md

This document elaborates on the cornerstones of the project. While the main Readme.md is for guidance of the visitor, this is to explain the

* objective,
* deliverables (when we have all these, the paper can be submitted),
* literature relevant to the study, with a short description,
* and relevant communication, e.g. from people who do not contribute to the repo.

### 02-06_Abc.md

Here, the manuscript can be drafted. In 03\_Results.md, figures are implemented. 02_Methods and 03\_Results can contain a supplementary section.

### environment.yml

Ideally, the complete project can be run from the command line in a defined environment like conda, that takes care of all dependencies. Any environment or containerization service is possible, as seen fit. If graphical programs are used, whose results need to be exactly reproducible, the graphical program is started from the command line in the respective environment or container. In case of conda usage you can define a basic environment.yml that can be used to created the projects conda environment with `sh create_conda_env.sh`. Updates to an existing conda environment can be performed with `sh update_conda_from_env.sh` and the further modified conda environment can then mirrored to an updated environment.yml file with `sh update_environment.sh` to enable reproducable analysis steps.

### Notes.md

This file serves as a pinboard for the project. The general description of the project is in 01_Background.md. Notes.md is for preserving thoughts and ideas for later. It should also contain a RAID log to write down all risks, assumptions, issues or dependencies that occurred during analysis that have the potential to make your work invalid or would need you to rerun steps. This is not a Kanban board. An external Kanban software can be useful in phases of intense collaboration.

## Notebooks

Every person participating in the project needs to document their work in a comprehensible manner. There are numerous functional ways to log your work. However, collaboration is easier if everybody adheres to the same structure. Every person keeps the files in their own Notebook_YourName folder.

### Notebook README.md

A table of contents listing your log files and describing in few words what has been done. If any output has been generated, the name of the analysis step directory is noted, e.g. 02_berryreplace. This is useful to `grep` all info on one step.

### Single step doc files

Every logically connected analysis step gets one document. You can start by copying the notebookTemplate.md to your notebook folder. Some suggestions:

* Document dead ends and failures.
* The numbering in your notebook will not match the analysis step subdirectories. Keep track by mapping your logs to analysis steps in the README.md of your notebook folder.
* An analysis step might be documented in multiple docs, but no doc covers multiple analysis steps.

Every doc file contains  

* Title
* Date
* Computer and path where the step was executed.
* Objective and reason
* Commands run. Use wrapper scripts to simplify.
* Description of the result

Every sensible unit of work has its own doc file.  
If there are any concerns regarding risks, assumptions, issues or dependencies of this step, then note them down in Notes.md in the main directory.

## Analysis step subdirectories

An analysis step is the way from an input to an output that is an input to another step or a result. Examples include "Coexpression analysis", "Quality control", "Differential expression analysis". Each directory contains

* A wrapper bash-script which reruns the whole step.
* Scripts that are called in the wrapper.
* Directories local/, temp/ and/or results/.

temp/ is used for any files that can be safely deleted, but are kept to study the workings of the script or to minimize time for a rerun. Could also be called work/. local/ is for results that are used later on, but which should not be synchronized with github because of size or privacy. results/ is for results and files used later on that can be uploaded to github safely.

The goal of this structure is that you can run the wrapper script from inside the analysis folder and the step is completely reproduced.

### Dealing with intermediate files not in the repo

Files in local/ directories are not synchronized with the repo, so after downloading the repo, a step depending on a "local" file cannot simply rerun. Prior to publication, such data can either be deposited in a data storage accessible for collaborators or we put up with it and rerun all steps.

## Credit

 This is developed from the [repo skeleton of ISU genomics](https://github.com/ISUgenomics/Repo_skeleton) described in detail [here](https://bioinformaticsworkbook.org/projectManagement/Intro_projectManagement.html).
