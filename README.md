# README

This is a workflow for calling CNV in WES data, and is currently under active development.

It utilises Docker and WDL/Cromwell, with backend support for SGE and GoogleCloud.

## Docker

This workflow uses several Docker images, with the philosophy of one tool -> one image.
Dockerfiles are listed under each tool subdir in the Docker folder.

### Building Docker images

e.g. commands for building and testing Docker images


`docker build --no-cache -t joshmschmidt/wescnv:0.0.1 .`


`docker run --rm -ti joshmschmidt/wescnv:0.0.1  /bin/bash`



### Profiles

There are some basic profiles that can be selected, located in conf/

garvan.config provides a nice basic sge env template.

*TODO:*

dynamic resource, particulalry memory and TMPDIR on cluster.

`nextflow run wescnv.nf -resume -profile garvan --inputFile inputFiles.txt`
