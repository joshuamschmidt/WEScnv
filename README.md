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
