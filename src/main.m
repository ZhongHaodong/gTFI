%   Created by Haodong Zhong on 2024.11.25
%   Last modified by Haodong Zhong on 2024.12.02

%% GRE Data Preprocessing Example
% The basic variable names related to QSM follow the conventions of the MEDI_toolbox.
filepath{1} = '';
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,BW,DICOM_dir]=Read_DICOM(filepath{1});
iMag = sqrt(sum(abs(iField).^2,4));
matrix_size = size(iMag);
Mask_BET = BET(iMag, matrix_size, voxel_size);
[iFreq_raw, N_std] = Fit_ppm_complex_TE(iField, TE);
iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
R2s = arlo(TE, abs(iField));

%% Generate the Mesh

%% Example of gTFI using the FMM form of L
params_gTFI = [];
params_gTFI.delta_TE = delta_TE;
params_gTFI.CF = CF;
params_gTFI.voxel_size = voxel_size;
params_gTFI.mask = mask_gTFI;
params_gTFI.iMag = iMag;
params_gTFI.f = iFreq;
params_gTFI.N_std = N_std;
params_gTFI.B0_dir = B0_dir;
params_gTFI.lambda = 1250;
params_gTFI.R2s = R2s;
params_gTFI.Mask_CSF = Mask_CSF; 
%params_gTFI.flag_CSF=false;
[gTFI,RDF_gTFI] = gTFI_FMM(params_gTFI,mesh);
