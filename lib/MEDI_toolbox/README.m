% This is an implementation of the Morphology Enabled Dipole Inversion (MEDI)
% method for reconstructing a Quantitative Susceptibility Map from MR data.
% The code is not fully optimized and is given for educational purpose.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use this tool box, add MEDI_toolbox to your MATLAB Path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE DATASETS %%%%%%%%%%%%%%%%%%%%%%
% Example datasets can be found MEDI_data 
% 01_Numerical_phantom contains the simulation in Neuroimage 2012;59(3):2560-8.
% 02_Wienieff_Liu contains a numerical brain
% 03_Invivo_GE contains a human brain dataset acquired from a GE scanner
% 04_Invivo_Siemens contains a human brain dataset acquired from a Siemens scanner
%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE DATASETS %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%
%
% high dimensional variable
%
% iField - 4, or 5 dimensional complex MRI dataset. 
%          the 4th dimension is echo
%          the 5th dimension is channel
%
% 3D variables
%
% Mask - binary mask denoting the region of interest
% iMag - magnitude image, square root of squares of all echoes
% iFreq_raw - the raw field map, which may contain wrapping, 
%             unit in rad/echo
% N_std - estimated noise standard deviation on iFreq_raw
% iFreq - the unwrapped field map, aka total field
%         unit in rad/echo
% RDF - Relative Difference Field, aka local field
%       unit in rad/echo
% QSM - Quantitative Susceptibility Map, 
%       unit in parts per million, aka ppm
%
% vectors
%
% B0_dir - unit vector representing direction of B0 field 
% matrix_size - sizes ([x y z]) of the imaging volume
% voxel_size - size of the voxel
%              unit in mm
% TE - echo time, unit in sec
%
% scalars
%
% delta_TE - echo spacing, unit in sec
% CF - center frequency, unit in Hz
% B0_strength - magnetic field strength, unit in Tesla
%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% USEFUL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('functions\','functions\_LBV', 'functions\_spurs_gc');


% IMPORTANT: please copy <Your MATLAB installation path>\toolbox\images\iptformats\private\dicomparse.mexw64 to the 'functions' directory of the toolbox (confirm overwrite when asked).
% The following lines of code should take care of this for you!
% Use accelerated (for Siemens and GE only) reading of DICOMs
try 
    [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_xxx_DICOM('DICOM_dir');
    % xxx=GE or Siemens or Philips.
catch
    copyfile([matlabroot '\toolbox\images\iptformats\private\dicomparse.mexw64'], 'functions\dicomparse.mexw64')
end
try
    [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_xxx_DICOM('DICOM_dir');
catch
    [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_xxx_DICOM_old('DICOM_dir');
end
    
    
% In case of Bruker real/imag data folders:
[iField, CF, B0_dir, Affine3D, TE, delta_TE, matrix_size, voxel_size] = Read_Bruker_DICOM('real_dir','imag_dir');
% Other formats may be supported in the future

% Remove echo-to-echo phase inconsistencies in readout phase corected
% complex data
% NOTE: initial testing suggests that use of this function DOES NOT
% introduce negative effects further down the pipeline if normal data
% (i.e., without phase correction) is processed with this function.
[iField] = iField_correction(iField,voxel_size);

% Estimate the frequency offset in each of the voxel using a complex
% fitting (even echo spacing)
[iFreq_raw N_std] = Fit_ppm_complex(iField);

% Estimate the frequency offset in each of the voxel using a complex
% fitting (uneven echo spacing)
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField,TE);


% Spatial phase unwrapping (region-growing)
 iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);

% Spatial phase unwrapping (graph-cut based)
 iFreq = unwrapping_gc(iFreq_raw,iMag,voxel_size);
%%%% Simultaneous Phase Unwrapping and Removal of Chemical Shift (SPURS) Using Graph Cuts: Application in Quantitative Susceptibility Mapping
%%%% IEEE TME 20015;34(2):531-540

% if large fringe lines persists, try 
 iFreq = unwrapLaplacian((iFreq_raw, matrix_size, voxel_size);

 %Prepare mask based on magnitude thresholding
 Mask = genMask(iField,voxel_size);
 
% Background field removal using Projection onto Dipole Fields
 RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);
%%%% NMR Biomed 2011;24(9):1129-36.
%%%% MRM 2010;63(1):194-206

% Background field removal using Laplacian Boundary Value
 RDF = LBV(iFreq,Mask,matrix_size,voxel_size,0.005);
%%%% NMR Biomed 2014;27(3):312-319

% Before running MEDI, variables need to be saved
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir;

 QSM = MEDI_L1('lambda',1000);
% morphology enabled dipole inversion

% write QSM as DICOMs
 write_QSM_dir(QSM,'DICOM_dir','QSM_DICOM')

 Visu3D( QSM, 'dimension',voxel_size);
% Visualize a 3D matrix

write_QSM_dir(QSM, 'DICOM_dir', 'QSM_DICOM');
% Save results in DICOM format

% This is version 0.5a of the code. 
% We tried to make sure that all the functions that are needed are present
% However, if you are missing some, please contact qsmreconstruction@gmail.com.

% Enjoy!



