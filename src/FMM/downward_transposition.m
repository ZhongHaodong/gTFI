function [u, du] = downward_transposition(c_nm,p,params_FMM)

%% Initialization
OT = params_FMM.OcTree;
R_gq_f_MOM = params_FMM.R_gq_f_MOM;
dRdn_gq_f_MOM = params_FMM.dRdn_gq_f_MOM;
Interaction_List = params_FMM.Interaction_List;
S_conj = params_FMM.S_conj;
R_L2L = params_FMM.R_L2L;
binAssignment = params_FMM.binAssignment;

leaf_depth = max(OT.BinDepths);
d_nm = zeros([OT.BinCount,p+1,2*p+1]);
%% downword
cell_Level2 = find(OT.BinDepths==2);
for i = 1:length(cell_Level2)
    Interaction_i = Interaction_List(cell_Level2(i),:);
    Interaction_i = Interaction_i(~isnan(Interaction_i));
    lI = length(Interaction_i);
    if ~isempty(Interaction_i)
        for n = 0:p
            for m = -n:n
                S = S_conj(cell_Level2(i),1:lI,n+1:n+1+p,m+p+1:m+3*p+1);
                S = reshape(S,[lI,p+1,2*p+1]);
                d_nm(cell_Level2(i),n+1,m+p+1) = sum((-1)^n.*S.*c_nm(Interaction_i,:,:),'all');
            end
        end
    end
end

R_L2L_expans = padarray(R_L2L,[0,0,0,p],0);
if leaf_depth>3
    for cell_depth = 3:leaf_depth-1
        cell_childLevel = find(OT.BinDepths==(cell_depth+1));
        cell_ThisLevel = unique(OT.BinParents(cell_childLevel));
        c_map = zeros([length(cell_ThisLevel),189,p+1,2*p+1],'gpuArray');
        S_map = gpuArray(S_conj(cell_ThisLevel,:,:,:));
        for i = 1:length(cell_ThisLevel)
            Interaction_i = Interaction_List(cell_ThisLevel(i),:);
            Interaction_i = Interaction_i(~isnan(Interaction_i));
            if ~isempty(Interaction_i)
                c_map(i,1:length(Interaction_i),:,:) = c_nm(Interaction_i,:,:);
            end
        end
        for n = 0:p
            for m = -n:n
                d_M2L = sum((-1)^n.*S_map(:,:,n+1:n+1+p,m+p+1:m+3*p+1).*c_map(:,:,:,:),[2,3,4]);
                d_L2L = zeros([length(cell_ThisLevel),1]);
                for n2=n:p
                    R = R_L2L_expans(cell_depth,binAssignment(cell_ThisLevel),n2-n+1,2*p+1-m-n2:2*p+1-m+n2);
                    R = reshape(R,[length(cell_ThisLevel),2*n2+1]);
                    d = d_nm(OT.BinParents(cell_ThisLevel),n2+1,p+1-n2:p+1+n2);
                    d = reshape(d,[length(cell_ThisLevel),2*n2+1]);
                    d_L2L = d_L2L+sum(R.*d,2);
                end
                d_nm(cell_ThisLevel,n+1,m+p+1) = d_M2L+d_L2L;
            end
        end
    end
end
%% leaf
Batch = 4000;
cell_leafLevel = find(OT.BinDepths==leaf_depth);

cell_leafLevel = cell_leafLevel(~isnan(cell_leafLevel));
Batch_Num = ceil(length(cell_leafLevel)/Batch);
for j = 1:Batch_Num
    if j<Batch_Num
        c_map = zeros([Batch,189,p+1,2*p+1],'gpuArray');
        S_map = gpuArray(S_conj(cell_leafLevel(1+(j-1)*Batch:j*Batch),:,:,:));
        for i = 1:Batch
            Interaction = Interaction_List(cell_leafLevel(i+(j-1)*Batch),:);
            Interaction = Interaction(~isnan(Interaction));
            if ~isempty(Interaction)
                c_map(i,1:length(Interaction),:,:) = c_nm(Interaction,:,:);
            end
        end
        for n = 0:p
            for m = -n:n
                d_M2L = sum((-1)^n.*S_map(:,:,n+1:n+1+p,m+p+1:m+3*p+1).*c_map(:,:,:,:),[2,3,4]);
                d_L2L = zeros([Batch,1]);
                for n2=n:p
                    R = R_L2L_expans(leaf_depth,binAssignment(cell_leafLevel(1+(j-1)*Batch:j*Batch)),n2-n+1,2*p+1-m-n2:2*p+1-m+n2);
                    R = reshape(R,[Batch,2*n2+1]);
                    d = d_nm(OT.BinParents(cell_leafLevel(1+(j-1)*Batch:j*Batch)),n2+1,p+1-n2:p+1+n2);
                    d = reshape(d,[Batch,2*n2+1]);
                    d_L2L = d_L2L+sum(R.*d,2);
                end
                d_nm(cell_leafLevel(1+(j-1)*Batch:j*Batch),n+1,m+p+1) = d_M2L+d_L2L;
            end
        end

    else
        if mod(length(cell_leafLevel),Batch) ~= 0
            Mod = mod(length(cell_leafLevel),Batch);
        else
            Mod = Batch;
        end
        c_map = zeros([Mod,189,p+1,2*p+1],'gpuArray');
        S_map = gpuArray(S_conj(cell_leafLevel(1+(j-1)*Batch:end),:,:,:));
        for i = 1:Mod
            Interaction = Interaction_List(cell_leafLevel(i+(j-1)*Batch),:);
            Interaction = Interaction(~isnan(Interaction));
            if ~isempty(Interaction)
                c_map(i,1:length(Interaction),:,:) = c_nm(Interaction,:,:);
            end
        end

        for n = 0:p
            for m = -n:n
                d_M2L = sum((-1)^n.*S_map(:,:,n+1:n+1+p,m+p+1:m+3*p+1).*c_map(:,:,:,:),[2,3,4]);
                d_L2L = zeros([Mod,1]);
                for n2=n:p
                    R = R_L2L_expans(leaf_depth,binAssignment(cell_leafLevel(1+(j-1)*Batch:end)),n2-n+1,2*p+1-m-n2:2*p+1-m+n2);
                    R = reshape(R,[Mod,2*n2+1]);
                    d = d_nm(OT.BinParents(cell_leafLevel(1+(j-1)*Batch:end)),n2+1,p+1-n2:p+1+n2);
                    d = reshape(d,[Mod,2*n2+1]);
                    d_L2L = d_L2L+sum(R.*d,2);
                end
                d_nm(cell_leafLevel(1+(j-1)*Batch:end),n+1,m+p+1) = d_M2L+d_L2L;
            end
        end
    end
end
%% 分发
d_map = d_nm(OT.PointBins,:,:);
u = sum(R_gq_f_MOM.*d_map,[2,3]);
du = sum(dRdn_gq_f_MOM.*d_map,[2,3]);
end