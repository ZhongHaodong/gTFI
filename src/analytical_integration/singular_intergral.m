function I = singular_intergral(v1,v2,v3)
%   When using the code, please cite 
%   Qinlong Ren, Cho Lik Chan,
%   Engineering Analysis with Boundary Elements,
%   Volume 53,
%   Pages 1-8,
%   ISSN 0955-7997,
%   https://doi.org/10.1016/j.enganabound.2014.11.018.

barycentric = (v1+v2+v3)./3;
v12 = v1-barycentric;
v22 = v2-barycentric;
v32 = v3-barycentric;
eqs = 1e-12;
x=v12/sqrt(sum(v12.^2));
z=cross(x,v22);
y0 = cross(z,x);
y = y0./sqrt(sum(y0.^2));
v13 = zeros(2);
v23 = zeros(2);
v33 = zeros(2);

v13(1) = dot(v12,x);
v13(2) = dot(v12,y);
v23(1) = dot(v22,x);
v23(2) = dot(v22,y);
v33(1) = dot(v32,x);
v33(2) = dot(v32,y);

detA = (v13(1)-v33(1))*(v23(2)-v33(2))-(v13(2)-v33(2))*(v23(1)-v33(1));
if detA == 0
    detA = detA+eqs;
end
m11 = (v23(2)-v33(2))/detA;
m21 = -1*(v13(2)-v33(2))/detA;
m12 = -1*(v23(1)-v33(1))/detA;
m22 = (v13(1)-v33(1))/detA;
a=(m22^2+m21^2)/(m11*m22-m21*m12)^2;
b=(m12^2+m11^2)/(m11*m22-m21*m12)^2;
c=-2*(m22*m12+m21*m11)/(m11*m22-m21*m12)^2;
FF = 1/3*(log((4*a+2*b+2*sqrt(4*a+b-2*c)*sqrt(a+b-c)-3*c)/(-2*a-4*b+2*sqrt(a+4*b-2*c)*sqrt(a+b-c)+3*c))/sqrt(a+b-c)...
    +log((4*a+2*sqrt(a)*sqrt(4*a+b-2*c)-c)/(-2*a-c+2*sqrt(a)*sqrt(a+b+c)))/sqrt(a)...
    +log(((4*b+2*sqrt(b)*sqrt(a+4*b-2*c)-c)*(2*b+c+2*sqrt(b)*sqrt(a+b+c)))/((b+sqrt(b)*sqrt(4*a+b-2*c)-c)*(-1*b+sqrt(b)*sqrt(4*a+b-2*c)+c)))/sqrt(b));
J = (v13(1)*v23(2)+v23(1)*v33(2)+v33(1)*v13(2))-(v23(1)*v13(2)+v33(1)*v23(2)+v13(1)*v33(2));
I = 1/(4*pi)*FF*J;
end