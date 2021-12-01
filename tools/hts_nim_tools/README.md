# hts_nim_tools

This image is based on hts_nim_tools from Pederson: []


`docker build -t joshmschmidt/hts_nim_tools:0.2.0 .`

`docker push joshmschmidt/hts_nim_tools:0.2.0`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/hts_nim_tools:0.2.0  docker://joshmschmidt/hts_nim_tools:0.2.0 `

`docker run -u $(id -u):$(id -g) --rm -v $(pwd):/data joshmschmidt/hts_nim_tools:0.2.0 /bin/bash -c \
"hts_nim_tools count-reads \
--fasta Homo_sapiens_assembly38.fasta \
--mapq 20 \
gencode.v38.AUTO-ALL-EXONS-1TPG-MERGED-50bp_MS-GC.GC5-DF-SD.bed \
C6DAUANXX_3_ACG048_4.cram | sort -k1,1 -k2,2n | gzip > ACG048_4.cpt.txt.gz"`
