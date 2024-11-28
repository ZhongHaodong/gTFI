%校正奇偶回波线形相位差
%在图像域对相位差3D线形拟合
function iField_correction = iField_correction_bipolar(iField,voxel_size,mask)
iField_correction = iField;
matrix_size = size(squeeze(iField(:,:,:,1)));
mag = sqrt(sum(abs(iField).^2,4));
pha =2*angle(iField(:,:,:,2)) - (angle(iField(:,:,:,1)) + angle(iField(:,:,:,3)));

%Quick and dirty unwrapping
pha = unwrapPhase(mag,pha,matrix_size);
%pha = unwrapLaplacian(pha, matrix_size, voxel_size);

%Filter and mask pha
kernel_filter2 = ones(3)*(1/9);
pha_filter2 = zeros(matrix_size);
for i=1:matrix_size(3)
    pha_filter2(:,:,i) = filter2(kernel_filter2,pha(:,:,i));
end
mask_filter2 = abs(pha-pha_filter2) < 0.2;
Mask = logical(mask.*mask_filter2);
pha = pha(Mask);

%Create location matrix
X_Loc = zeros(matrix_size);
Y_Loc = zeros(matrix_size);
Z_Loc = zeros(matrix_size);
X_Loc = repmat((1:matrix_size(1))'*voxel_size(1),[1,matrix_size(2),matrix_size(3)]);
Y_Loc = repmat((1:matrix_size(2))*voxel_size(2),[matrix_size(1),1,matrix_size(3)]);
zloc = zeros(1,1,matrix_size(3));
zloc(1,1,:) = (1:matrix_size(3))*voxel_size(3);
Z_Loc = repmat(zloc,[matrix_size(1),matrix_size(2),1]);

X_Loc_filter2 = X_Loc(Mask);
Y_Loc_filter2 = Y_Loc(Mask);
Z_Loc_filter2 = Z_Loc(Mask);
Unit = ones(size(X_Loc_filter2));
A = [Unit X_Loc_filter2 Y_Loc_filter2 Z_Loc_filter2];

y = A\pha;

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
         temp = temp.*exp(-1i*0.25*pha0);
    else
         temp = temp.*exp(1i*0.25*pha0);
    end
    iField_correction(:,:,:,i) = reshape(temp,matrix_size);
end
end