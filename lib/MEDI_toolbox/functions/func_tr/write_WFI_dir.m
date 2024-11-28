%   Last modified by Alexey Dimov on 2016.05.13

function y = write_WFI_dir(WFI,dicomDir,saveDir,type)
warning( 'off', 'all' );

if nargin<3
    if ~exist('WFI_DICOM','dir')
        mkdir('WFI_DICOM')
    end
    saveDir = 'WFI_DICOM';
end

% if ~exist('qsm_DICOM','dir')
%     mkdir('qsm_DICOM')
% end
WFI = permute(WFI,[2,1,3]); %change row/column order due to differences in representations in DICOM and Matlab
WFI = int16(WFI*1000);

filelist = dir(dicomDir);
i=1;
while i<=length(filelist)
    if filelist(i).isdir==1
        filelist = filelist([1:i-1 i+1:end]);   % eliminate folders
    else
        i=i+1;
    end
end

sliceNum = size(WFI,3);
sliceIdx = 0;
imIdx = 1; 
UID = dicomuid;
saveDir = [saveDir '/' type];
mkdir(saveDir);

while sliceIdx < sliceNum
    fid = fopen([dicomDir '/' filelist(imIdx).name]);
    if fid>0
        sliceIdx=sliceIdx + 1;
        fclose(fid);
        info = dicominfo([dicomDir '/' filelist(imIdx).name]);
        %--
        info.SeriesDescription = type;
        info.SeriesInstanceUID = UID;
        info.SeriesNumber = 99;
        info.EchoNumber = 1;
        info.EchoTime = 0.0;
        info.InstanceNumber = sliceIdx;
        %--
        dicomwrite(WFI(:,:,sliceIdx),[saveDir '/' num2str(sliceIdx) '.dcm'], info);
    end
    imIdx = imIdx + 1;
end
addpath(saveDir);
warning( 'on', 'all' )