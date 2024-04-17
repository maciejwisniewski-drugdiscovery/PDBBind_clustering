#!/bin/bash
#SBATCH --job-name=distance_matrix
#SBATCH --nodes=1
#SBATCH --cpus-per-task=38
#SBATCH --mem=200G
#SBATCH --time=24:00:00

input_fasta=/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Cluster/PDBBind_proteins_sequences.fasta
distant_matrix_output=/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Cluster/clustalo_PDBBind_protein_sequences_dsitance_matrix.txt


cd /mnt/evafs/groups/sfglab/mwisniewski/anaconda3/bin
. activate
source activate PDBBind_clustering
cd /mnt/evafs/groups/sfglab/mwisniewski/ingenix/PDBBind_clustering

clustalo -i $input_fasta --distmat-out=$distant_matrix_output --full --force --threads=38
