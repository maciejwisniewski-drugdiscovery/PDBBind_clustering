#!/bin/bash
#SBATCH --job-name=cdhit_clusters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=24:00:00

input_fasta=/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta
output_clusters=/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/cdhit_PDBBind_protein_sequences_clusters.clstr


cd /mnt/evafs/groups/sfglab/mwisniewski/anaconda3/bin
. activate
source activate PDBBind_clustering
cd /mnt/evafs/groups/sfglab/mwisniewski/ingenix/PDBBind_clustering

cd-hit -i $input_fasta -o $output_clusters -c 0.9 -n 5 -d 0 -T 20
