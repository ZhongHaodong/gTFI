function params_FMM = Tree_parameters(params_FMM, mesh,p,voxel_size,gq_number)
OT = params_FMM.OcTree;
Mask = OT.Mask;
nf = size(mesh.faces,1);
TR = triangulation(mesh.faces,double(mesh.vertices));
nv = faceNormal(TR);
Interaction_List = params_FMM.Interaction_List;

cell_center = (((OT.BinIndex_Range(:,1:3)+ OT.BinIndex_Range(:,4:6))./2)...
    -(1+size(Mask))./2).*(voxel_size)';
cell_center = single(cell_center);
O = cell_center(OT.PointBins(:),:);
[Xq, W] = quadrature(mesh, gq_number); % Quadrature points and weights
%% upword
O_gq = repmat(O,[1,1,gq_number]);
O_gq = permute(O_gq,[3,1,2]);
O_gq = single(reshape(O_gq,[gq_number*nf,3]));

nv_gq = repmat(nv,[1,1,gq_number]);
nv_gq = permute(nv_gq,[3,1,2]);
nv_gq = single(reshape(nv_gq,[gq_number*nf,3]));

R_gq_f_MOM = zeros([nf,p+1,2*p+1],'single','gpuArray');
dRdn_gq_f_MOM = zeros([nf,p+1,2*p+1],'single','gpuArray');
for n = 0:p
    [R_nm_f,dRdn_nm_f] = SH_ExpansionR_dRdn(O_gq, Xq, nv_gq, n);
    R_gq_f_MOM(:,n+1,p+1-n:p+1+n) = squeeze(sum(reshape(W.*R_nm_f,[gq_number,nf,2*n+1]),1));
    dRdn_gq_f_MOM(:,n+1,p+1-n:p+1+n) = squeeze(sum(reshape(W.*dRdn_nm_f,[gq_number,nf,2*n+1]),1));
end
R_M2M = zeros([OT.BinCount,p+1,2*p+1],'single','gpuArray');
for n = 0:p
    R_M2M(2:end,n+1,p+1-n:p+1+n) = SH_ExpansionR(cell_center(OT.BinParents(2:end),:),cell_center(2:end,:), n);
end
%% downword
R_L2L_map = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
R_L2L = zeros(max(OT.BinDepths),8,p+1,2*p+1);
binAssignment = zeros([OT.BinCount,1],'int8');
for i = 2:OT.BinCount
    gtMask = bsxfun(@gt, cell_center(i,:),cell_center(OT.BinParents(i),:));
    [~,binAssign] = max(all(bsxfun(@eq,gtMask,R_L2L_map),2),[],3);
    binAssignment(i) = binAssign;
end

for depth = 2:(max(OT.BinDepths)-1)
    cell_ThisLevel = find(OT.BinDepths==depth);
    for i = 1:length(cell_ThisLevel)
        child = find(OT.BinParents==cell_ThisLevel(i));
        if length(child)==8
            for j = 1:8
                binAssign = binAssignment(child(j));
                for n = 0:p
                    R_L2L(depth+1,binAssign,n+1,p+1-n:p+1+n) = SH_ExpansionR(cell_center(cell_ThisLevel(i),:),cell_center(child(j),:),n);
                end
            end
        end
    end
end

leaf_size = size(Mask)./2^max(OT.BinDepths);
[Y,X,Z]=meshgrid(-(leaf_size(2)-1)/2:(leaf_size(2)-1)/2,...
    -(leaf_size(1)-1)/2:(leaf_size(1)-1)/2,...
    -(leaf_size(3)-1)/2:(leaf_size(3)-1)/2); %cell_map
X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);
X_list = reshape(X,[numel(X),1]);
Y_list = reshape(Y,[numel(X),1]);
Z_list = reshape(Z,[numel(X),1]);
R_final0 = zeros(numel(X),p+1,2*p+1,'single');

for n = 0:p
    R_final0(:,n+1,p+1-n:p+1+n) = SH_ExpansionR(zeros(numel(X),3),[X_list,Y_list,Z_list],n);
end
R_final = reshape(R_final0,[leaf_size(1),leaf_size(2),leaf_size(3),n+1,2*p+1]);
S = zeros([OT.BinCount,189,2*p+1,4*p+1],'single');

for i = 2:OT.BinCount
    Interaction_i = Interaction_List(i,:);
    Interaction_i = Interaction_i(~isnan(Interaction_i));
    lI = length(Interaction_i);
    cc = cell_center(Interaction_i,:);
    if ~isempty(Interaction_i)
        pn = repmat(cell_center(i,:),[length(Interaction_i),1]);
        Si = zeros(189,2*p+1,4*p+1);
        for n = 0:2*p
            Si(1:lI,n+1,2*p+1-n:2*p+1+n) = SH_ExpansionS(cc, pn, n);
        end
        S(i,:,:,:) = Si;
    end
end
S_conj = conj(S);

%% Assign variables to fields in params_FMM
params_FMM.R_gq_f_MOM = R_gq_f_MOM;
params_FMM.dRdn_gq_f_MOM = dRdn_gq_f_MOM;
params_FMM.cell_center = cell_center;
params_FMM.R_M2M = R_M2M;
params_FMM.S_conj = S_conj;
params_FMM.R_L2L = R_L2L;
params_FMM.R_final = R_final;
params_FMM.binAssignment = binAssignment;
fprintf('Tree_parameters completed.\n')
end

