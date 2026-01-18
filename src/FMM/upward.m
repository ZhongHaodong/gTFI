function c_nm = upward(fs,p,params_FMM)
% Upward : function that computes the moments on all levels starting from
%          the leaves up to the second level. The moments are approximated
%          up to p terms. Leaves have a direct computation whereas parent
%          cells' moments are computed by an M2M translation from children.
%
% INPUTS
% ------------------------------------------------- regarding info and data
% - Xq      : solution at previous step
% - p   : number of terms in the multipole expansion

% OUTPUTS
% - c_nm      : multipole expansion moments of every cell (matrix nexp x ncells)

%% Initialization
OT = params_FMM.OcTree;
R_gq_f_MOM = params_FMM.R_gq_f_MOM;
dRdn_gq_f_MOM = params_FMM.dRdn_gq_f_MOM;
R_M2M = params_FMM.R_M2M;
D2N = params_FMM.D2N;
c_nm = zeros(OT.BinCount,p+1,2*p+1);
u_n = D2N*fs.*gather(R_gq_f_MOM)-fs.*gather(dRdn_gq_f_MOM);
leaf_depth = max(OT.BinDepths);

%% leaf_moment

cell_ThisLevel = find(OT.BinDepths==leaf_depth);
for cell = 1:length(cell_ThisLevel)
    c_nm(cell_ThisLevel(cell),:,:) = sum(u_n(OT.PointBins==cell_ThisLevel(cell),:,:),1);
end

%% M2M
for cell_depth = (leaf_depth-1):-1:2
    c_nm_expansion = gpuArray(padarray(c_nm,[0,0,p],0));
    cell_ThisLevel = find(OT.BinDepths==cell_depth);
    cell_childLevel = find(OT.BinDepths==(cell_depth+1));
    BinParents = OT.BinParents;
    M2M = zeros(length(cell_childLevel),p+1,2*p+1);

    for n = 0:p
        for m = -n:n
            M2M(:,n+1,m+p+1) = squeeze(sum(R_M2M(cell_childLevel,1:n+1,:).*c_nm_expansion(cell_childLevel,n+1:-1:1,((m+2*p+1)+p:-1:(m+2*p+1)-p)),[2,3]));
        end
    end

    for i = 1:length(cell_ThisLevel)

        cell_child = find(BinParents==cell_ThisLevel(i));
        if isempty(cell_child)== 0
            c_nm(cell_ThisLevel(i),:,:) = squeeze(sum(M2M(ismember(cell_childLevel,cell_child),:,:),1));
        end
       
    end

end

end