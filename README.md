Convolutional bayesian models for anatomical landmarking on multi-dimensional shapesGaussian Process based Landmarking (MICCAI2020)
By Yonghui Fan, -7/30/2022

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

Please cite:
Fan Y, Wang Y. Convolutional Bayesian Models for Anatomical Landmarking on Multi-Dimensional Shapes. Med Image Comput Comput Assist Interv. 2020;12264:786-796. doi: 10.1007/978-3-030-59719-1_76. Epub 2020 Sep 29. PMID: 34291235; PMCID: PMC8291336.
