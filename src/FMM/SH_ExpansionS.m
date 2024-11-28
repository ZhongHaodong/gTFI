function S = SH_ExpansionS(o, x, n)
% |ox|>|oy|
    %ox = gpuArray(x-o);
    ox = x-o;
    ox_sph = cartesian_to_spherical(ox);
    %P_ox_pos = gpuArray(permute(legendre(n,cos(ox_sph(:,2))),[2,1]));
    P_ox_pos = permute(legendre(n,ox_sph(:,2)),[2,1]);
    m = 1:n;
    P_ox_neg = (-1).^m.*factorial(n-m)./factorial(n+m).*P_ox_pos(:,2:end);
    m = -n:n;
    %S = factorial(n-m).*[P_ox_neg,P_ox_pos].*exp(1i.*m.*ox_sph(:,3)).*ox_sph(:,1).^(-(n+1));
    S = factorial(n-m).*[P_ox_neg(:,end:-1:1),P_ox_pos].*exp(1i.*m.*ox_sph(:,3)).*(1./ox_sph(:,1).^(n+1));
    
end
function spherical_coor = cartesian_to_spherical(xyz)
%transform cartesian coordinates (x, y, z) to spherical
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

xy = x.^2 + y.^2;

r = sqrt(xy + z.^2);
%theta = acos(z./r);
cos_theta = z./r;
phi = atan2(y, x);
spherical_coor = [r, cos_theta, phi];
end