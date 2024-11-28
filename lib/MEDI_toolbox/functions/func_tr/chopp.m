%a  = iField(:,:,:,:);
function iField_out = chopp(iField)
 iField_out =iField;
for i =1:60
    if mod(i,2)==0
        iField_out(:,:,i,:) =   iField(:,:,i,:) .*exp(1i*pi);
    end
end
end