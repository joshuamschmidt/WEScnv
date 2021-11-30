# mosdepth

This image is based on mosedpth from Pederson and Quinlan: []

Built form static binary:

`docker build -t joshmschmidt/mosdepth:0.3.2 .`

`docker push joshmschmidt/mosdepth:0.3.2`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/mosdepth:0.3.2 docker://joshmschmidt/mosdepth:0.3.2`

`docker run -u $(id -u):$(id -g) --rm -v $(pwd):/data joshmschmidt/mosdepth:0.3.2 /bin/bash -c \
"mosdepth --fasta Homo_sapiens_assembly38.fasta \
--by gencode.v38.AUTO-ALL-EXONS-1TPG-MERGED-50bp_MS-GC.GC5-DF-SD.bed \
--no-per-base \
--mapq 20 \
--threads 1 \
ACG048_4 \
C6DAUANXX_3_ACG048_4.cram"`
