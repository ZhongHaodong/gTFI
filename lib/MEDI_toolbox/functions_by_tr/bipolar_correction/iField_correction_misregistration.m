function iField_interpolation = iField_correction_misregistration(iField,voxel_size,Mask,BW,f0)
%% Default readout direction =row
%% iFreq = 2pi * f0 * delta_TE
[sx sy sz necho] = size(iField);
% ²åÖµ±¶ÊýL
L =100; 
FOVx = sy *voxel_size(1);
iField_interpolation = zeros([sx sy sz necho]);
xLoc = zeros(sy,1);
Linex = zeros(sy,1); 
InterpoLinex = zeros(L*sy,1);

% remove  iFreq constant offset  
   a =f0(Mask);
   a =a(a~=0);
   m = mean(a);
   f0 = f0-m;
   
for i=1:sx
    for j=1:sz
        for k=1:necho
            xinitial = int32((1:sy)*L) ;
            xshift =  (-1)^k * int32(f0(i,:,j) / BW * L);    %unit pixel
           %xshift =1;
            xLoc = xinitial + xshift;
            xLoc(xLoc<=0 | xLoc>= sy * L) = xinitial(xLoc<=0 | xLoc>= sy * L);
            Linex = iField(i,:,j,k);
            InterpoLinex = interpft(Linex , sy * L);
            try
             iField_interpolation(i,:,j,k) = InterpoLinex(xLoc);
            catch
                xLoc
            end
        end
    end
end
end