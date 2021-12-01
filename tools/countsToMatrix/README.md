
`conda create --name jossch_countsToMatrix python=3.10 pandas numpy`

`conda activate jossch_countsToMatrix`

`pip freeze > requirements.txt`

`docker build -t joshmschmidt/counts_to_matrix:0.0.1 .`

`docker push joshmschmidt/counts_to_matrix:0.0.1`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/counts_to_matrix:0.0.1  docker://joshmschmidt/counts_to_matrix:0.0.1 `
