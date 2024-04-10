# Old Matlab Codes

The codes provided in this repo are **not maintained** anymore. However, from time to time, I need to share some of them for a project, hence here they are. 

Most likely useful : 
* *Anisotropic Diffusion using Lattice Basis Reduction.*
The subdirectory AD_LBR of this repo constructs the sparse matrix of this method (Fehrenbach, Mirebeau, JMIV, 2013). This is a sparse non-negative scheme for divergence form diffusion, with positive definite but possibly strongly anisotropic diffusion tensors. Alternatively, see [this notebook](https://github.com/Mirebeau/AdaptiveGridDiscretizations/blob/main/Notebooks_Div/AnisotropicDiffusion.ipynb) for a Python implementation with some explanations and illustrations. 


The other files from this repo are unlikely to be useful without some input from me. 

Regarding *Anisotropic eikonal equation solvers*, please use my more recent and maintained repositories, which contain a few Matlab interfaces  : 
[AdaptiveGridDiscretizations](https://github.com/Mirebeau/AdaptiveGridDiscretizations/tree/main/Notebooks_FMM/Matlab), and [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching/tree/master/Interfaces/MatlabHFM).