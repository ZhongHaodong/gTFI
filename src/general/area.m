function s = area(v1,v2,v3)
% Calculate the areas or triangles
v12 = v1-v2;
v13 = v1-v3;
C = cross(v12,v13,2);
s = 1/2.*(sum(C.^2,2)).^0.5;
end