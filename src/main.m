%   Created by Haodong Zhong on 2024.11.25
%   Last modified by Haodong Zhong on 2024.12.02

%% GRE Data Preprocessing Example
% The basic variable names related to QSM follow the conventions of the MEDI_toolbox.
filename{1} = 'E:\MRI_DATA\BEM_TFI\health_final\DATA_TOTAL\1_BEM_ZHD20240715M\QSM_tra_mono_ 1mm_6echo';
[iField0,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,BW,DICOM_dir]=Read_DICOM_qsmliver([filename{1},'\data']);
iField = iField0;
iMag = sqrt(sum(abs(iField(:,:,:,4:end)).^2,4));
matrix_size = size(iMag);
Mask = BET(iMag, matrix_size, voxel_size);
[iFreq_raw, N_std] = Fit_ppm_complex_TE(iField, TE);
iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
R2s = arlo(TE, abs(iField));

%% Generate the Mesh

%% Example of gTFI using the FMM form of L
params_gTFI = [];
params_gTFI.delta_TE = delta_TE;
params_gTFI.CF = CF;
params_gTFI.voxel_size = voxel_size;
params_gTFI.mask = mask;
params_gTFI.iMag = iMag;
params_gTFI.f = iFreq;
params_gTFI.N_std = N_std;
params_gTFI.B0_dir = B0_dir;
params_gTFI.lambda = 1000;
params_gTFI.R2s = R2s;
params_gTFI.Mask_CSF = logical(extract_CSF(R2s, mask, voxel_size));

gTFI = gTFI_FMM(params_gTFI,mesh);

export_folder = [filename{1},'\results_gTFI'];
mkdir(export_folder);
nii = make_nii(fliplr(gTFI),voxel_size,[],16);
save_nii(nii,[export_folder,'\gTFI.nii']);
