function fs_FMM = LH_FMM(fb,p,params_FMM,near_field)
c_nm = upward_transposition(fb,p,params_FMM);
[fs_far, dfs_far] = downward_transposition(c_nm,p,params_FMM);
fs_FMM = real(params_FMM.D2NH*(gather(fs_far/(4*pi))+near_field.T'*double(fb(:))) ...
    -(gather(dfs_far/(4*pi))+near_field.Q'*double(fb(:))));
end