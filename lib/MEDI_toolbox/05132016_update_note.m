% This is an updated version of MEDI toolbox.
% 
Important changes:

1. New function to read GE DICOMs - significant speedup is achieved. Read_GE_DICOM now properly handles datasets acquired using ZIP option on scanner. New function uses the same interface as the previous verision.

Note: since the new function uses additional files which might potentially be a source of compatibility issues between different versions of Matlab, old function is still included in the toolbox as "Read_GE_DICOM_old".

2. Fit_ppm_complex allows estimation of field map from single echo datasets. Single echo data can now be processed with the same reconstruction pipeline.

3. Fixed bug in write_QSM_dir (used to generate DICOM series for reconstructed QSM) which previously caused improper slice location.

Please contact qsmreconstruction@gmail.com with your questions and
suggestion

Regards,
Cornell MRI lab