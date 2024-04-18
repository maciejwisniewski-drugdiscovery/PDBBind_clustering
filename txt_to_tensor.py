import torch
import numpy as np
from numpy import genfromtxt


distant_matrix_txt='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/clustalo_PDBBind_protein_sequences_distance_matrix.txt'
tensor_distance='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/clustalo_PDBBind_protein_sequences_distance_matrix.pt'
tensor_distance_gpu='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/gpu_clustalo_PDBBind_protein_sequences_distance_matrix.pt'

numpy_matrix = np.genfromtxt(distant_matrix_txt, skip_header=1)

numpy_matrix = numpy_matrix[:,1:]

torch_tensor = torch.from_numpy(numpy_matrix)
torch_tensor_gpu = torch_tensor.to(device='cuda')

torch.save(torch_tensor, tensor_distance)
torch.save(torch_tensor_gpu, tensor_distance_gpu)