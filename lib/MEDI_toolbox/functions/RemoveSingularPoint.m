% Remove singular point due to failure in Fat/Wat seperation
%   [iFreq_SingularPointRemoved Mask_EdgeErased] = RemoveSingularPoint(iFreq, Mask,
%   matrix_size)

% 
%   output
%   iFreq_SingularPointRemoved -
%
%   input
%   iFreq - 
%   Mask - 
%   matrix_size - dimension of the field of view
%   
%   MATLAB  created by Jianqi Li(jqli@phy.ecnu.edu.cn) on 2015.03.05
%   


function [iFreq_SingularPointRemoved, Mask_EdgeErased] = RemoveSingularPoint(iFreq, Mask, matrix_size)

disp('removing the singular point in frequency map ....');

Mask_EdgeErased = Mask;
iFreq_SingularPointRemoved = iFreq;
IsSingularPoint = zeros(matrix_size);
fPhaseMean = zeros(matrix_size);

%The first Slice
k = 1;
for j = 2:matrix_size(2)-1;
    for i= 2:matrix_size(1)-1;
         Mask_EdgeErased(i,j,k) = 0;
    end
end

%The last Slice
k = matrix_size(3);
for j = 2:matrix_size(2)-1;
    for i= 2:matrix_size(1)-1;
         Mask_EdgeErased(i,j,k) = 0;
    end
end

for k = 2:matrix_size(3)-1;
    for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if( (Mask(i,j,k) + Mask(i-1,j,k) + Mask(i+1,j,k) + Mask(i,j-1,k) + Mask(i,j+1,k) + Mask(i,j,k-1) + Mask(i,j,k+1)) == 7)
                Mask_EdgeErased(i,j,k) = 1;
            else
                Mask_EdgeErased(i,j,k) = 0;
            end
        end
    end
end

%The first Slice
k = 1;
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k+1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>2.0)
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
end

%The last Slice
k = matrix_size(3);
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k-1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>2.0)
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
end

     
for k = 2:matrix_size(3)-1;
    for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k)+ IsSingularPoint(i,j,k-1)+ IsSingularPoint(i,j,k+1);
                 fPhaseMean(i,j,k) = fPhaseSum / (6 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>2.0)
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
     end
end

%The first Slice
k = 1;
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k+1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>1.0)
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
end

%The last Slice
k = matrix_size(3);
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k-1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>1.0)
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
end

for k = 2:matrix_size(3)-1;
    for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k)+ IsSingularPoint(i,j,k-1)+ IsSingularPoint(i,j,k+1);
                 %奇异点不能与平均值计算
                 if((Mask(i,j,k)==1)&&(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>1.0))
                      IsSingularPoint(i,j,k)=1;               
                 end
            end
         end
     end
end

%The first Slice
k = 1;
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k+1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>0.5)
                      IsSingularPoint(i,j,k)=1; 
                      iFreq_SingularPointRemoved(i,j,k) = fPhaseMean(i,j,k);
                 end
            end
         end
end

%The last Slice
k = matrix_size(3);
for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k) + IsSingularPoint(i,j,k-1);
                 fPhaseMean(i,j,k) = fPhaseSum / (5 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>0.5)
                      IsSingularPoint(i,j,k)=1;  
                      iFreq_SingularPointRemoved(i,j,k) = fPhaseMean(i,j,k);
                 end
            end
         end
end

for k = 2:matrix_size(3)-1;
    for j = 2:matrix_size(2)-1;
        for i= 2:matrix_size(1)-1;
            if(Mask(i,j,k)==1)
                fPhaseSum =    iFreq(i-1,j,k)*(1-IsSingularPoint(i-1,j,k)) + iFreq(i+1,j,k)*(1-IsSingularPoint(i+1,j,k))+ iFreq(i,j-1,k)*(1-IsSingularPoint(i,j-1,k)) + iFreq(i,j+1,k)*(1-IsSingularPoint(i,j+1,k)) + iFreq(i,j,k-1)*(1-IsSingularPoint(i,j,k-1)) + iFreq(i,j,k+1)*(1-IsSingularPoint(i,j,k+1));
                nTotlaSingularPoint = IsSingularPoint(i-1,j,k)  + IsSingularPoint(i+1,j,k) + IsSingularPoint(i,j-1,k)+ IsSingularPoint(i,j+1,k)+ IsSingularPoint(i,j,k-1)+ IsSingularPoint(i,j,k+1);
                 fPhaseMean(i,j,k) = fPhaseSum / (6 - nTotlaSingularPoint); 
                 %奇异点不能与平均值计算
                 if((Mask(i,j,k)==1)&&(abs(iFreq(i,j,k) - fPhaseMean(i,j,k))>0.5))
                      IsSingularPoint(i,j,k)=1; 
                      iFreq_SingularPointRemoved(i,j,k) = fPhaseMean(i,j,k);
                 end
            end
         end
     end
end

end
