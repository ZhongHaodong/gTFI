function fb_FMM = L_FMM(fs,p,params_FMM,near_field)
OT = params_FMM.OcTree;

c_nm =upward(fs,p,params_FMM);
fb_far = downward(c_nm,p,params_FMM);
%fb_FMM =real(reshape(near_field.T*(params_FMM.D2N*fs) ...
%    -near_field.Q*fs,size(OT.Mask)));
fb_FMM =real(fb_far/(4*pi) + reshape(near_field.T*(params_FMM.D2N*double(fs)) ...
    -near_field.Q*double(fs),size(OT.Mask)));
end