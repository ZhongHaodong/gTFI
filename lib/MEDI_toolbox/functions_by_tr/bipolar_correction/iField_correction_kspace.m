%Zero and First order Phase error correction between odd and even echos for brain in kspace 
%Need Search_rho_plus.m
%Created by Rui Tong
function iField_out= iField_correction_kspace(iField,voxelSize,mis_direct,mask)

% mis_direct =
% 1  misalignment occurs in column direction .
% 2 misalignment occurs in row direction .          1 default.

if nargin < 4
    mask = genMask(iField,voxelSize);
end
mask = logical(mask);

if nargin < 3
    mis_direct =1;
end

iField_out = iField;
row = size(iField,1);
col = size(iField,2);
slice = size(iField,3);
NumEcho = size(iField,4);

%标准化，用元胞数组存数据
Data_image = cell(slice ,NumEcho);
Data_kspace = cell(slice ,NumEcho);
for j=1:slice
    for i=1:NumEcho
        Data_image{j,i} =iField(:,:,j,i);                                % i 表示回波，j表示层数 
        if (mis_direct == 2)                                        %若一阶相位差出现在图像行方向，则置换行列
             Data_image{j,i} = Data_image{j,i}.';
        end
        Data_kspace{j,i} = fftshift(fft2(Data_image{j,i}));    %对图像执行fft2，得到k空间数据
    end
      %E(:,:,j) = angle(Data_image{j,2}) -0.5 *(angle(Data_image{j,1}) + angle(Data_image{j,3}));     %展示偶回波与奇回波的相位差
end

% if (mis_direct == 2)                                     
%     mask = mask';
% end
% for i =1:6
%     F{i} = D{16,i}(:,84);
%     figure(i)
%     plot(abs(F{i}));
% end
middle_slice = int16(slice/2);    % 选取中心层
row = size(Data_kspace{1,1},1);
col = size(Data_kspace{1,1},2);

%% first order correction with correlation
    phase = @(f,m,alpha,beta)              exp(1i* (alpha *m + beta * f));          %%附加相位函数

    %options = optimset('PlotFcns',@saplotf,'MaxIter',200);
    options = optimset('MaxIter',200);
     u_left = -2;
     u_right = 2;
     F = repmat([0:row-1]',1,col);
     M = ones(row,col);
    [ibeta, rho1,exitflag,output] = simulannealbnd(@(x)Search_rho_plus(x,middle_slice,Data_kspace,F), 0,u_left, u_right,options)  %  模拟退火
    
 %% 利用ibeta执行一阶相位校正
    %beta =0;
    beta = ibeta*2*pi/ size(F,1)
    alpha = 0*pi;
    for j=1:slice
        for i=1:NumEcho            
            if mod(i,2) ==0
                %C{j,i} = C{j,i}.*phase(F,M,alpha,beta);
                Data_image{j,i} =  ifft2(fftshift(ifft(fft(Data_kspace{j,i}).*phase(F,M,alpha,beta))));            
%             else 
%                  Data_image{j,i} =  ifft2(fftshift(ifft(fft(Data_kspace{j,i}).*phase(F,M,alpha,-beta))));
            end
            if (mis_direct == 2)
                Data_image{j,i} = Data_image{j,i}.';
            end
            %Data_image{j,i} = Data_image{j,i}.*mask(:,:,j);
            %H{j,i} = fftshift(fft2(Data_image{j,i}));
        end
    end 
    
%     %% 在图像域进行零阶相位校正
%     for j=1:size(iField,3)
%         G(:,:,j) = Data_image{j,1};
%         %E(:,:,j) = angle(B(:,:,j,4)) -0.5 *(angle(B(:,:,j,5)) + angle(B(:,:,j,3))); 
%         E(:,:,j) = angle(exp(2i*(angle(Data_image{j,4}) -0.5 *(angle(Data_image{j,5}) + angle(Data_image{j,3}))))); 
%         %E(:,:,j) = angle(Data_image{j,2}) - angle(sqrt(Data_image{j,1}.*Data_image{j,3}));
%     end
%     
%     E_mask = E(mask);
%     phi0 = angle(sum(exp(1i*E_mask(:))))/2;
%      for i=1:NumEcho
%          if mod(i,2) ==0
%              for j = 1:slice 
%             	Data_image{j,i} = Data_image{j,i} * exp(-1i*phi0);
%              end
%          end
%      end
    
     for j=1:slice
        for i=1:NumEcho
           iField_out(:,:,j,i) =Data_image{j,i}; 
        end
     end

end
