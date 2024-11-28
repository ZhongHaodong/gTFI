% This is an updated version of MEDI toolbox.
% 
Important changes:
1. Forward difference with Neuman boundary conditions and backward
difference with Dirichlet boundary conditions are now used to compute
discrete gradient and discrete divergence (functions fgrad and bdiv).
This approach suppresses "checkerboard pattern" artifact in final QSM.
2. GUI version is now available

New functions:
1. Read_Bruker_DICOM and Read_Bruker_raw functions are now available to
load Bruker data.
2. Fit_ppm_complex_TE function to estimate frequency offset for datasets
with uneven echo spacing
3. unwrap_gc - graph-cut based phase unwrapping
4. LBV - background field removal by solving Laplacian boundary value
problem
5. write_QSM_dir - outputs reconstructed QSM as a set of DICOM images

Please contact qsmreconstruction@gmail.com with your questions and
suggestion

Regards,
Cornell MRI lab