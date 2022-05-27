###Create a custome reference to be used for cellranger count.sh
##inputs: Fasta file and gtf file
###outputs: genome directory with custom name

#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=70:00:00


cd /home/des65576/4913-345456114/

ml CellRanger/4.0.0
cellranger mkref --genome=stickleback_may25 --fasta=stickleback_hicGenome.fa --genes=10XscRNAd.gtf
mkref (END)

