%函数功能：校正奇偶回波线形相位差
%方法：ErrorMap是由iField生成的奇偶回波相位差图，在图像域对ErrorMap3D线形拟合，所得参数来校正iField

function iField_correction = iField_correction_ErrorMap_3DLinear(iField,ErrorMap,voxel_size,mask)
iField_correction = iField;
matrix_size = size(squeeze(iField(:,:,:,1)));
mag = sqrt(sum(abs(iField).^2,4));

%Quick and dirty unwrapping
ErrorMap1 =unwrapping_gc(ErrorMap,mag,voxel_size);%Haodong
%ErrorMap1 = unwrapPhase(mag,ErrorMap1,matrix_size);
%pha = unwrapLaplacian(pha, matrix_size, voxel_size);

%Filter and mask pha
ErrorMap_filtered = zeros(matrix_size);
for i=1:matrix_size(3)
   ErrorMap_filtered(:,:,i) = medfilt2(ErrorMap1(:,:,i),[5 5]);
   %ErrorMap_filtered(:,:,i) = filter2( fspecial('average',3),ErrorMap1(:,:,i));
end
mask_filtered = abs(ErrorMap1- ErrorMap_filtered) > 0.2;
m1 = SMV(mask_filtered,matrix_size, voxel_size, 3)>0.001;
Mask = logical(mask.*(~m1));
% Mask1 = double(m1);
% Mask2 = double(Mask);
% ErrorMap1(~Mask) = 0;
ErrorMap1 = ErrorMap_filtered(Mask);

X_Loc = zeros(matrix_size);
Y_Loc = zeros(matrix_size);
Z_Loc = zeros(matrix_size);
X_Loc = repmat((1:matrix_size(1))',[1,matrix_size(2),matrix_size(3)]);
Y_Loc = repmat((1:matrix_size(2)),[matrix_size(1),1,matrix_size(3)]);
zloc = zeros(1,1,matrix_size(3));
zloc(1,1,:) = (1:matrix_size(3));
Z_Loc = repmat(zloc,[matrix_size(1),matrix_size(2),1]);

X_Loc_filtered = X_Loc(Mask);
Y_Loc_filtered = Y_Loc(Mask);
Z_Loc_filtered = Z_Loc(Mask);
Unit = ones(size(X_Loc_filtered));
A = [Unit X_Loc_filtered Y_Loc_filtered Z_Loc_filtered];

y = A\ErrorMap1;
fprintf('3DLinear fit results:\n  intercept: %f\n  slope1:%f\n  slope2:%f\n  slope3:%f\n',y(1),y(2),y(3),y(4));

X_Loc = X_Loc(:);
Y_Loc = Y_Loc(:);
Z_Loc = Z_Loc(:);
Unit = ones(size(X_Loc));
A = [Unit X_Loc Y_Loc Z_Loc];
pha0=A*y;

for i=1:size(iField,4)
    temp = iField(:,:,:,i);
    temp = temp(:);
    if mod(i,2) == 0        
         temp = temp.*exp(-1i*pha0);
    else
         temp = temp.*exp(1i*pha0);
    end
    iField_correction(:,:,:,i) = reshape(temp,matrix_size);
end
end