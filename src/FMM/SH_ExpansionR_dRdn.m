function [R,dRdn] = SH_ExpansionR_dRdn(o, y, nv, n)% 1 OK
% |ox|>|oy|
if n==0
    R = ones(size(y,1),1);
    dRdn = zeros(size(y,1),1);
else
    oy = gpuArray(y-o);
    oy_sph = cartesian_to_spherical(oy);
    P=legendre(n,cos(oy_sph(:,2)));
    P_oy_pos = gpuArray(permute(P,[2,1]));
    dP_oy_pos = gpuArray(permute(legendre_derivative(n,P,cos(oy_sph(:,2))),[2,1]));
    m = 1:n;
    P_oy_neg = (-1).^m.*factorial(n-m)./factorial(n+m).*P_oy_pos(:,2:end);
    dP_oy_neg = (-1).^m.*factorial(n-m)./factorial(n+m).*dP_oy_pos(:,2:end);
    m = -n:n;

    P_oy = [P_oy_neg(:,end:-1:1),P_oy_pos];
    dP_oy = [dP_oy_neg(:,end:-1:1),dP_oy_pos];
    R = (1./factorial(n+m)).*P_oy.*exp(1i*m.*oy_sph(:,3)).*oy_sph(:,1).^n;
    
    dRdr = (1./factorial(n+m)).*P_oy.*exp(1i*m.*oy_sph(:,3)).*(n*oy_sph(:,1).^(n-1));
    dRdTheta = (-1./factorial(n+m)).*exp(1i*m.*oy_sph(:,3)).*dP_oy.*sin(oy_sph(:,2)).*oy_sph(:,1).^(n-1);% without -sin(Sita)
    dRdPhi = 1i*m.*(1./factorial(n+m)).*P_oy.*exp(1i*m.*oy_sph(:,3)).*oy_sph(:,1).^(n-1);

    drdnv = dot(oy./oy_sph(:,1),nv,2);
    dThetadnv = cos(oy_sph(:,3)).*cos(oy_sph(:,2)).*nv(:,1)...
        + sin(oy_sph(:,3)).*cos(oy_sph(:,2)).*nv(:,2)...
        - sin(oy_sph(:,2)).*nv(:,3);
    dPhidnv = dot([-1*oy(:,2),oy(:,1),zeros(size(oy(:,1)))]./sqrt(oy(:,1).^2+oy(:,2).^2),nv,2);
    dRdn = dRdr.*drdnv + dRdTheta.*dThetadnv + dRdPhi./sin(oy_sph(:,2)).*dPhidnv;

end
end
function spherical_coor = cartesian_to_spherical(xyz)
%transform cartesian coordinates (x, y, z) to spherical
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

xy = x.^2 + y.^2;

r = sqrt(xy + z.^2);
theta = acos(z./r);
phi = atan2(y, x);
spherical_coor = [r, theta, phi];
end

