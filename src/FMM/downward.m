function u = downward(c_nm,p,params_FMM)

%% Initialization
OT = params_FMM.OcTree;
Interaction_List = params_FMM.Interaction_List;
S_conj = params_FMM.S_conj;
R_L2L = params_FMM.R_L2L;
R_final = params_FMM.R_final;
binAssignment = params_FMM.binAssignment;

leaf_depth = max(OT.BinDepths);
d_nm = zeros([OT.BinCount,p+1,2*p+1]);
%% downword
cell_childLevel = find(OT.BinDepths==3);
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
        d_nm(cell_ThisLevel,n+1,m+p+1) = d_M2L;
    end
end

R_L2L_expans = padarray(R_L2L,[0,0,0,p],0);
if leaf_depth>3
    for cell_depth = 3:1:leaf_depth-1
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
Batch = 10000;
cell_leafLevel = find(OT.BinDepths==leaf_depth);
for i = 1:length(cell_leafLevel)
    Index_Range = OT.BinIndex_Range(cell_leafLevel(i),:);
    cell_mask = OT.Mask(Index_Range(1):Index_Range(4),Index_Range(2):Index_Range(5),Index_Range(3):Index_Range(6));
    if all(cell_mask==0)
        cell_leafLevel(i) = NaN;
    end
end
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
d_final = d_nm(cell_leafLevel,:,:);
u = zeros(size(OT.Mask));
R_final = gpuArray(R_final);
for i =1:length(cell_leafLevel)
    Index_Range = OT.BinIndex_Range(cell_leafLevel(i),:);
    cell_mask = OT.Mask(Index_Range(1):Index_Range(4),Index_Range(2):Index_Range(5),Index_Range(3):Index_Range(6));
    if all(cell_mask==0)
        continue;
    end
    cell_mask = gpuArray(logical(cell_mask));
    d_map = repmat(d_final(i,:,:),[numel(cell_mask),1,1]);
    d_map = gpuArray(reshape(d_map,[size(cell_mask),p+1,2*p+1]));
    u(Index_Range(1):Index_Range(4),Index_Range(2):Index_Range(5),Index_Range(3):Index_Range(6)) =...
        sum(cell_mask.*R_final.*d_map,[4,5]);
end
end