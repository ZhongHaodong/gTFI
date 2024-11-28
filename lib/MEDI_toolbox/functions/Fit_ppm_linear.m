% Projection onto Dipole Fields (PDF)
%   [p1, dp1, relres, p0]=Fit_ppm_complex(M)
%    
%   output
%   p1 - field map, may need further unwrapping
%   dp1 - a priori error estimate
%   relres - relative residual
%   p0 - initial phase
%   Y - unwrapped iField
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   The coil combination method is similar to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Last modified by Alexey Dimov on 2016.05.12


function [p1, dp1, p0, Y]=Fit_ppm_linear(S,TE)

%Modification to handle one echo datasets - assuming zero phase at TE = 0;
%- AD, 05/12/2016
if size(S,4) == 1
    S = cat(4,abs(S),S);
end

if size(S,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    S = sum(S.*conj( repmat(S(:,:,:,1,:),[1 1 1 size(S,4) 1])),5);  
    S = sqrt(abs(S)).*exp(1i*angle(S));
end


S= conj(S);
s0=size(S);
L_s0=length(s0);
nechos=size(S,L_s0);

S=reshape(S,[prod(s0(1:L_s0-1)),s0(L_s0)]);

Y=angle(S);
c=((Y(:,2)-Y(:,1)));
[m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
c(ind==1)=c(ind==1)-2*pi;
c(ind==3)=c(ind==3)+2*pi;

for n=1:nechos-1
    cd=((Y(:,n+1)-Y(:,n)))-c;
    Y(cd<-pi,(n+1):end)=Y(cd<-pi,n+1:end)+2*pi;
    Y(cd>pi,(n+1):end)=Y(cd>pi,n+1:end)-2*pi;
end
A = [1  0 ;1 1; 1 2 ;1 3;1 4;1 5;1 6;1 7];
ip = A(1:nechos,:)\Y(:,1:nechos)';
p0 = ip(1,:)';
p1 = ip(2,:)';

% dp1 = p1;
 dp1 = ones(size(p1));
% iMag = sqrt(sum(abs(S).^2,2));
% dp1 = 1./iMag;

p1=reshape(p1,s0(1:L_s0-1));
dp1=reshape(dp1,s0(1:L_s0-1));
p0=reshape(p0,s0(1:L_s0-1));
Y=reshape(Y,[s0(1:L_s0-1) 8]);
    

