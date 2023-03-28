# Developer's guide

This document is a brief overview of how to use this code base. Some understanding of Nextflow and how it implements DSL2 is assumed. 

## Basic concept

This pipeline base is organized in the following way:

* `main.nf` - entry point into the pipeline, imports the core workflow from workflow/<pipeline>.nf
* `workflow/<pipeline.nf>` - the actual core logic of the pipeline; imports sub-workflows from subworkflow/<sub>.nf
* `subworkflow/<sub>.nf` - a self-contained processing chain that is part of the larger workflow (e.g. read alignment and dedup in a WGS calling workflow)
* `modules/<module>.nf` - A command line tool/call that can be imported into a (sub)workflow. 

## Groovy libraries

This pipeline imports a few functions into the nextflow files from lib/ - mostly to keep the actual pipeline code a bit cleaner/more readable. For example, 
the `--help` command line option can be found in lib/WorkflowMain.groovy. Likewise, you could use this approach to do some basic validation of your inputs etc. 

## Github workflows

Github supports the automatic execution of specific tasks on code branches, such as the automatic and building pushing of Docker containers. To add github workflows to your
repository, rename the folder `dot_github` to `.github` and adjust the files therein accordingly (name of pipeline, docker repo etc).

### Docker containers

In order to automatically push Docker containers, you must add your docker username and API token as secrets to your repository (DOCKERHUB_USERNAME and DOCKERHUB_TOKEN). Secrets 
can be created under Settings/Secrets and Variables/Actions. Of course, you also need to have an account on Dockerhub. 







