####Aligns fastqs to genome directory, counts cells and creates gene expression matrix. Output can be exported into loupe, seurat etc.
##Fastqs must be in direcotry with path listed
###Requires original --id name
### Ensure gtf is formatted correctly and contains transcript and gene ID for each line

#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=30
#SBATCH --nodes=1
#SBATCH --mem=18gb
#SBATCH --time=168:00:00


cd /scratch/des65576/singlecellmay/

ml CellRanger/4.0.0

cellranger count --id=stick_Testis_may25v5 \
--fastqs=/home/des65576/4913-345456114/fastqs/ \
--sample=4913_Stickleback_testis_1,4913_Stickleback_testis_2,4913_Stickleback_testis_3,4913_Stickleback_testis_4 \
--transcriptome=/home/des65576/4913-345456114/stickleback_v5c \
--chemistry=SC3Pv3
