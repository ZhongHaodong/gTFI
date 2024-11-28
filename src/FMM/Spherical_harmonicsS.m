function YS = Spherical_harmonicsS(o, y, n)
% |ox|>|oy|
if n==0
    YS = ones(size(y,1),1);
else
    %oy = gpuArray(y-o);
    oy = y-o;
    oy_sph = cartesian_to_spherical(oy);
    P_oy_pos = permute(legendre(n,cos(oy_sph(:,2))),[2,1]);
    P_oy_neg = P_oy_pos(:,end:-1:2);
    m = -n:n;
    YS = (-1).^m.*sqrt(factorial(n-abs(m))./factorial(n+abs(m))).*[P_oy_neg,P_oy_pos].*exp(1i.*m.*oy_sph(:,3)).*oy_sph(:,1).^-(n+1);

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
