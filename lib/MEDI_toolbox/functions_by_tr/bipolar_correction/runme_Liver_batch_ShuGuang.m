%clear
import_folder = [];
tempFolder  = cell(0);
path = 'F:\MRIdata\TongRui\QSM_Liver\1ShuGuangHosp\910';
filelist = dir(path);
i=1;
    while(i<=length(filelist))
        if(filelist(i).isdir == 1)
            tempFolder = [tempFolder,filelist(i)];
        end
        i =i +1;
    end
    tempFolder = tempFolder(3:numel(tempFolder));
    numFolderPath = numel(tempFolder);
     filename = cell(numFolderPath,1);
   for i = 1:numFolderPath
       filename{i} = [path,'\',tempFolder{i}.name];
   end  

 %filename{1} ='F:\MRIdata\TongRui\QSM_Liver\bipolar\new\liuxiaolong_20180814'
 %filename{2}='F:\MRIdata\TongRui\QSM_Liver\1ShuGuangHosp\songnianbao';
 %filename{3}='F:\MRIdata\TongRui\QSM_Liver\曙光 医院\QSM HUANG\VIBE QSM SL2.5mm';
% filename{1}= 'F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM汪建';
  %filename{1}= 'F:\MRIdata\TongRui\QSM_Liver\1ShuGuangHosp\chentianbiao';
%   filename{2} ='F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM朱世杰';
  %filename{1} = 'F:\MRIdata\TongRui\QSM_Liver\zhaoyu_liver_20180502\vibe_q-dixon_tra_p4_bh_10echo_bip';
  % filename{1}='F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM Duhuifang';
   %filename{2}='F:\MRIdata\TongRui\QSM_Liver\zhaoyu_liver_20180502\vibe_q-dixon_tra_p4_bh_12echo_bip';
   %filename{4}= 'F:\MRIdata\TongRui\QSM_Liver\zhaoyu_liver_20180409\QSM_tra_p4_bip_deltaTE1.15';
   %filename{3}='F:\MRIdata\TongRui\QSM_Liver\zhaoyu_liver_20180409\vibe_QSM_tra_p4_mono_deltaTE1.79'
   %filename{1} = 'F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM 岑美英';
    %filename{2} = 'F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM 潘行建';
   % filename{1} = 'F:\MRIdata\TongRui\QSM_Liver\曙光医院\5.17 QSM 数据\5.17 QSM 数据\QSM 邓桂香';
    %filename{2} = 'F:\MRIdata\TongRui\QSM_Liver\曙光医院\QSM 张欢';
for i=1:4
    import_folder = [filename{i},'\mag_pha'];
   %import_folder = [filename{i},'\images_mag_phase_all_echoes']; 
format long;

[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,BW,DICOM_dir]=Read_DICOM_qsmliver(import_folder);
% TE = TE(1:6);
% iField = iField(:,:,:,1:6);
disp('Done Loading Dicom data');

% if i==3 || i==9
% SUBSAMPLE=2;
% if SUBSAMPLE == 2
%     iField = iField(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:,:);
%     voxel_size(1) = voxel_size(1)/0.5;
%     voxel_size(2) = voxel_size(2)/0.5;
% end
% end

iField = permute(iField,[2 1 3 4]);
voxel_size = voxel_size([2 1 3]);

 Mask = genMask_ECNU_2mm_3T_Liver(iField, voxel_size);
%Mask= genMask_Liver(iField, voxel_size,0.5);
mask = double(Mask);
matrix_size = size(Mask);

if (~exist('noise_level','var'))
    noise_level = calfieldnoise(squeeze(iField), Mask);
end

iField = iField/noise_level;

iMag = sqrt(sum(abs(iField).^2,4));

dfat = -3.5e-6*CF;

% 
    [wwater wfat wfreq R2s ErrorMap] = fit_IDEAL_R2_bip(conj(iField), TE, CF);
    iField_corr =  iField_correction_ErrorMap_3DLinear(conj(iField),ErrorMap,voxel_size,Mask);  
   [water fat iFreq R2s unwph_uf unwph N_std] = spurs_gc_R2s(iField_corr,TE,CF,voxel_size);
%   [water fat iFreq R2s unwph_uf unwph N_std] = spurs_gc_R2s(conj(iField),TE,CF,voxel_size);

PDFF = abs(fat)./(abs(fat) + abs(water));
if mean(PDFF(Mask&(PDFF<1)&(PDFF>0)))>0.5
    PDFF = 1-PDFF;
    temp = fat;
    fat =water;
    water = temp;
    clear temp;
end
%  [water_LM,fat_LM,r2star_LM] =fitLM_CPP(iField_corr,TE,abs(water),abs(fat),R2s,Mask,CF);
% %[water_LM,fat_LM,r2star_LM] =fitLM_CPP(iField_corr,TE,abs(fat),abs(water),R2s,Mask,CF);
% PDFF_LM = abs(fat_LM)./(abs(fat_LM) + abs(water_LM));

%N_std=ones(matrix_size);
[iFreq_SingularPointRemoved] = RemoveSingularPoint(iFreq, Mask, matrix_size);

%[RDF shim] = PDF(-iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);
%[RDF shim] = PDF(iFreq_SingularPointRemoved, N_std, Mask,matrix_size,voxel_size, B0_dir);
[RDF shim] = PDF(-iFreq_SingularPointRemoved, N_std, Mask,matrix_size,voxel_size, B0_dir);

% %%%% run MEDI %%%%%
% morphology enabled dipole inversion

%% write QSM as DICOMs
%  write_QSM_dir(QSM,DICOM_dir);
%export_folder = [filename{i},'\results_R2s_lambda1000_withCSF\water_fat_iFreq'];
export_folder = [filename{i},'\results_R2s_lambda1000\water_fat_iFreq'];
if ~exist(export_folder,'dir')
        mkdir(export_folder)
end
saveDir = export_folder;
% 
%
 water_nii = make_nii(rot90(abs(water),3));
 save_nii(water_nii,[saveDir,'\water_nii.img']);
 
  fat_nii = make_nii(rot90(abs(fat),3));
 save_nii(fat_nii,[saveDir,'\fat_nii.img']);
 
iFreq_nii = make_nii(rot90(iFreq,3));
 save_nii(iFreq_nii,[saveDir,'\iFreq_nii.img']);

PDFF_nii = make_nii(rot90(PDFF,3),voxel_size,[],16);
save_nii(PDFF_nii,[saveDir,'\PDFF_nii.img']);

R2s((R2s>1500)|(R2s<(-1500)))=0;
R2s_nii = make_nii(rot90(R2s,3),voxel_size,[],16);
save_nii(R2s_nii,[saveDir,'\R2s_nii.img']);


Mask_CSF = PDFF>0.8;
Mask_CSF = Mask_CSF & Mask; 

% STEP 3: Dipole Inversion 
% save RDF.mat RDF iFreq  iMag N_std Mask matrix_size...
%      voxel_size delta_TE CF B0_dir Mask_CSF;
% QSM = MEDI_L1L0('lambda',1000, 'lambda_CSF',100','smv',5);

save RDF.mat RDF iFreq  iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir ;
QSM = MEDI_L1('lambda',1000,'smv',5);

 QSM_nii = make_nii(rot90(QSM,3),voxel_size,[],16);
 %save_nii(QSM_nii,[saveDir,'\QSM_withCSF.img']);
  save_nii(QSM_nii,[saveDir,'\QSM.img']);
 if 1
       saveDir2 = [filename{i},'\results_R2s_lambda1000\qsm_DICOM'];
        write_QSM_dir(permute(QSM,[2 1 3]),DICOM_dir,saveDir2);
 end
 
 if 0
        [water_LM,fat_LM,r2star_LM] =fitLM_CPP(iField_corr,TE,abs(water),abs(fat),R2s,Mask,CF);
        PDFF = abs(fat_LM)./(abs(fat_LM) + abs(water_LM));
        export_folder2 = [filename{i},'\results_R2s_lambda1000\LM'];
        if ~exist(export_folder2,'dir')
                mkdir(export_folder2)
        end
        saveDir = export_folder2;

        water_nii = make_nii(rot90(abs(water_LM),3));
         save_nii(water_nii,[saveDir,'\water_nii.img']);

        fat_nii = make_nii(rot90(abs(fat_LM),3));
        save_nii(fat_nii,[saveDir,'\fat_nii.img']);

        PDFF_nii = make_nii(rot90(PDFF,3));
        save_nii(PDFF_nii,[saveDir,'\PDFF_nii.img']);

        R2s_nii = make_nii(rot90(r2star_LM,3));
        save_nii(R2s_nii,[saveDir,'\R2s_nii.img']);
        
       saveDir2 = [filename{i},'\results_R2s_lambda1000\WFI_DICOM'];
       write_WFI_dir(permute(PDFF,[2 1 3]),DICOM_dir,saveDir2,'PDFF');
       write_WFI_dir(permute(water_nii.img,[2 1 3]),DICOM_dir,saveDir2,'water');
       write_WFI_dir(permute(fat_nii.img,[2 1 3]),DICOM_dir,saveDir2,'fat');
       write_WFI_dir(permute(R2s_nii.img,[2 1 3]), DICOM_dir,saveDir2,'R2s');
 end
end

  