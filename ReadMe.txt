Gaussian Process based Landmarking
By Yonghui Fan, -4/12/2021

landmark_<dataset name>.m is the code for landmarking of certain datasets.

surface mesh dataset name = 	McGill
								ModelNet
								MPIFaust
								mt1
								oasis
								radius
								Schelling
								SHREC
								teeth
								TOSCA

tetrahedron mesh dataset name =     ADNI2
									Autism

BNN: number of neighborhood, a hyperparameter to be set manually, usually from 20 - 500, it depends on your data.

Some codes are for visualization, users can comment them off.

My kernel function is implemented in /toolbox/@Mesh/GetGPLmk_Fan.m 