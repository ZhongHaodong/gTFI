function [Frequency phi0] = FrequencyEst(M,TE)
s0=size(M);
L_s0=length(s0);
nechos=size(M,L_s0);
numvox=prod(s0(1:L_s0-1));
Y = reshape(M,[numvox,nechos]);
A = ones(1,nechos);
A=[A;TE']';
ip = A\Y';
Frequency = ip(2,:)';
phi0 = ip(1,:)';
Frequency = reshape(Frequency,s0(1:L_s0-1))*(TE(2)-TE(1));
phi0= reshape(phi0,s0(1:L_s0-1))*(TE(2)-TE(1));
end