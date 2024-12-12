function gTFI = gTFI_FMM(params_gTFI,mesh)
%   Created by Haodong Zhong on 2024.11.26
%   Last modified by Haodong Zhong on 2024.11.28
%% 必要参数
delta_TE = params_gTFI.delta_TE;
CF = params_gTFI.CF;
voxel_size = params_gTFI.voxel_size;
mask = params_gTFI.mask;
iMag = params_gTFI.iMag;
f = params_gTFI.f;
N_std = params_gTFI.N_std;
B0_dir = params_gTFI.B0_dir;
lambda = params_gTFI.lambda;
%% 可选参数
R2s = getOrDefault(params_gTFI, 'R2s', []);
flag_CSF = getOrDefault(params_gTFI,'flag_CSF', true);
if flag_CSF
    Mask_CSF = params_gTFI.Mask_CSF;
    fprintf('CSF regularization used\n');
    lam_CSF = getOrDefault(params_gTFI,'lam_CSF', 10);
end
max_iter = getOrDefault(params_gTFI, 'max_iter', 6);
tol_norm_ratio = getOrDefault(params_gTFI, 'tol_norm_ratio', 0.05);
%% Default parameters that are generally not modified
master_mode = getOrDefault(params_gTFI, 'master_mode', false);
if master_mode
    fprintf('Master mode enabled\n');
end
p =  getParam_Master(params_gTFI,'multipoleOrder', 5, master_mode);

R2sThresh = getParam_Master(params_gTFI,'R2sThresh', 40, master_mode);
PB =  getParam_Master(params_gTFI,'PB', 30, master_mode);

cg_max_iter = getParam_Master(params_gTFI, 'cg_max_iter', 600, master_mode);
cgTFI_LogInterval = getParam_Master(params_gTFI, 'cgTFI_LogInterval', ceil(cg_max_iter/4), master_mode);
cg_tol = getParam_Master(params_gTFI, 'cg_tol', 0.05, master_mode);
cg_tol_init = getParam_Master(params_gTFI, 'cg_tol_init', 0.01, master_mode);
cg_max_init = getParam_Master(params_gTFI, 'cg_max_init', 30, master_mode);
fullTFIMode = getParam_Master(params_gTFI, 'fullTFIMode', false, master_mode);
merit = getParam_Master(params_gTFI, 'merit', true, master_mode);
%% matrix size 的预处理
matrix_size = size(f);
old_matrix_size=matrix_size;
flag_matrix_size=0;
flag_size1=matrix_size(1)/16;
flag_size2=matrix_size(2)/16;
flag_size3=matrix_size(3)/16;
if (mod(matrix_size(1),16)==0)&(mod(matrix_size(2),16)==0)&(mod(matrix_size(3),16)==0)
    if  ((mod(flag_size1,2)==0)||(mod(flag_size2,2)==0)||(mod(flag_size3,2)==0))
        flag_matrix_size=1;
    end
end
if ~flag_matrix_size
    new_matrix_size=matrix_size;
    if mod(matrix_size(1),16)~=0
        new_matrix_size(1)=matrix_size(1)+(16-mod(matrix_size(1),16));
    end
    if mod(matrix_size(2),16)~=0
        new_matrix_size(2)=matrix_size(2)+(16-mod(matrix_size(2),16));
    end
    if mod(matrix_size(3),16)~=0
        new_matrix_size(3)=matrix_size(3)+(16-mod(matrix_size(3),16));
    end
    if  ((mod(flag_size1,2)==1)&(mod(flag_size2,2)==1)&(mod(flag_size3,2)==1))
        new_matrix_size(1)=new_matrix_size(1)+16;
    end
    matrix_size=new_matrix_size;
    size_change=ceil((matrix_size-old_matrix_size)/2);
    mask = padarray(mask,size_change);
    iMag = padarray(iMag,size_change);
    f = padarray(f,size_change);
    R2s = padarray(R2s,size_change);
end
%% FMM预计算
params_FMM = [];
% Calculate OCtree
nv = faceNormal(triangulation(mesh.faces, double(mesh.vertices)));
centroid = (mesh.vertices(mesh.faces(:,1),:) + ...
    mesh.vertices(mesh.faces(:,2),:) + ...
    mesh.vertices(mesh.faces(:,3),:)) ./ 3; % Triangle centroids
params_FMM.OcTree = OcTree_matrix(mask,centroid,voxel_size,4);
%计算仅依赖OCtree的参数
[params_FMM.Neighbor_list,params_FMM.Interaction_List] = neighbor_searcher(params_FMM);
params_FMM = Tree_parameters(params_FMM, mesh,p,voxel_size,7);

% Compute D2N
params_FMM.D2N = computeD2N(mesh, params_FMM, 7);
params_FMM.D2NH = (params_FMM.D2N)';
% Compute near-field integrals
near_field = [];
[near_field.T, near_field.Q] = cal_near_field(mesh,params_FMM,7,voxel_size);

%% Preconditioner
MaskP = SMV(mask,matrix_size, voxel_size, 1)>0.001;
P = ones(matrix_size);
if isfield(params_gTFI, 'R2s') && ~isempty(params_gTFI.R2s)
    P(R2s > R2sThresh) = PB;
else
    fprintf('Preconditioned process did not proceed as expected due to missing R2s.\n');
end

P(~mask) = PB;
P(~MaskP) = 1;
%%
c2g = @(x) gpuArray(x);
g2c = @(x) gather(x);
%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
data_weighting_mode = 1;
gradient_weighting_mode = 1;
grad = @fgrad;
div = @bdiv;
%%
N_std = N_std.*mask;
tempn = double(N_std);
D=dipole_kernel(matrix_size, voxel_size, B0_dir);
m = dataterm_mask(data_weighting_mode, tempn, mask);
wG = gradient_mask(gradient_weighting_mode, iMag, mask, grad, voxel_size);

res_norm_ratio_x = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);

e=0.0000001; %a very small number to avoid /0
e = e.*(2*pi*delta_TE*single(CF)/single(1e6))^2;
e = e.*P.*P;
%% CSF regularization
if flag_CSF
    LT_reg = @(x) Mask_CSF.*(x - mean(x(Mask_CSF)));
end
%% 待拟合量及拟合函数的初始化
x_length = matrix_size(1)*matrix_size(2)*matrix_size(3);
fb_length = size(nv,1);
total_length = x_length + fb_length;

x = zeros([total_length ,1]);
zero_fb = zeros([fb_length,1]);

Col2Map = @(x) reshape(x(1:x_length),matrix_size);
Map2Col =  @(x) reshape(x,[x_length,1]);

DconvL = @(dx) real(ifftn(D.*fftn(P.*Col2Map(dx)))...
    +L_FMM(dx(x_length+1:end),p,params_FMM,near_field));

DconvL_H = @(dx) [Map2Col(real(P.*ifftn(D.*fftn(dx))));...
    LH_FMM(dx,p,params_FMM,near_field)];

Dconv = @(dx) real(ifftn(D.*fftn(dx)));

L = @(x) L_FMM(x,p,params_FMM,near_field);

L_H = @(x) LH_FMM(x,p,params_FMM,near_field);

A_init= @(x) L_H(mask.*mask.*L(x));

%计算fs的初始值
x(x_length+1:end) = real(pcg(A_init,L_H(mask.*f),cg_tol_init,cg_max_init,[],[]));
fprintf('initial fb calculation completed.\n');
%% 高斯牛顿法
iter=0;
b0 = m.*exp(1i*f);
fprintf('Gauss-Newton iteration starts.\n');
while (res_norm_ratio_x>tol_norm_ratio)&&(iter<max_iter)
    tic
    iter=iter+1;
    if mod(iter, 2) == 1 || fullTFIMode
        Vr = 1./sqrt(abs(wG.*grad(real(P.*Col2Map(x)),voxel_size)).^2+e);
        w = m.*exp(1i*(ifftn(D.*fftn(P.*Col2Map(x)))+L(x(x_length+1:end))));
        reg = @(dx) [Map2Col(P.*div(wG.*(Vr.*(wG.*grad(real(P.*Col2Map(dx)),voxel_size))) ...
            ,voxel_size));zero_fb];

        if flag_CSF
            reg_CSF = @(dx) lam_CSF.*...
                [Map2Col(P.*LT_reg(LT_reg(real(P.*Col2Map(dx)))));zero_fb];
            reg = @(dx) reg(dx) + 2*reg_CSF(dx);
        end

        fidelity = @(dx) DconvL_H(conj(w).*w.*DconvL(dx));

        A =  @(dx) reg(dx) + 2*lambda*fidelity(dx);

        b = reg(x) + 2*lambda*DconvL_H( real(conj(w).*conj(1i).*(w-b0)) );

        if(iter==1)
            dx = real(cgsolve(A, -b, 0.01, cg_max_iter,cgTFI_LogInterval));
        else
            dx = real(cgsolve(A, -b, cg_tol, cg_max_iter,cgTFI_LogInterval));
        end
        %dx = real(gmres(A,-b,150,0.03,2));

        dx_print = [Map2Col(P.*Col2Map(dx));dx(x_length+1:end)];
        x_print = [Map2Col(P.*Col2Map(x));x(x_length+1:end)];
        res_norm_ratio = norm(dx_print(:))/norm(x_print(:));
        res_norm_ratio_x = norm(dx_print(1:x_length))/norm(x_print(1:x_length));
        res_norm_ratio_b = norm(dx_print(x_length:end))/norm(x_print(x_length:end));

        x = x + dx;

        wres=m.*exp(1i*(real(ifftn(D.*fftn(P.*Col2Map(x))))+...
            L(x(x_length+1:end)))) - b0;

        cost_data_history(iter) = norm(wres(:),2);
        cost=abs(wG.*grad(real(P.*Col2Map(x(1:x_length))),voxel_size));
        cost_reg_history(iter) = sum(cost(:));

        elapsed_time_sec = toc;
        minutes = floor(elapsed_time_sec / 60);
        seconds = floor(mod(elapsed_time_sec, 60));

        fprintf(['iter: %d;\n res_norm_ratio:%8.4f; res_norm_ratio_x:%8.4f; res_norm_ratio_b:%8.4f;\n' ...
            'cost_L2:%8.4f; cost:%8.4f;\n elapsed_time: %d:%02d (min:sec).\n'], ...
            iter, res_norm_ratio, res_norm_ratio_x, res_norm_ratio_b, ...
            cost_data_history(iter), cost_reg_history(iter), minutes, seconds);


    else %加速chi收敛
        x_inner = Col2Map(x);

        bi = m.*exp(1i*(f-L(x(x_length+1:end))));
        Vr = c2g(1./sqrt(abs(wG.*grad(real(P.*x_inner),voxel_size)).^2+e));
        w = c2g(m.*exp(1i*Dconv(P.*x_inner)));

        reg = @(dx) P.*div(wG.*(Vr.*(wG.*grad(real(P.*dx),voxel_size))),voxel_size);
        if flag_CSF
            reg_CSF = @(dx) lam_CSF.*P.*LT_reg(LT_reg(real(P.*dx)));
            reg = @(dx) reg(dx) + reg_CSF(dx);
        end

        fidelity = @(dx) P.*Dconv(conj(w).*w.*Dconv(P.*dx) );

        A =  @(dx) reg(dx) + 2*lambda*fidelity(dx);

        b = reg(x_inner) + 2*lambda*P.*Dconv( real(conj(w).*conj(1i).*(w-c2g(bi))) );

        dx = g2c(real(cgsolve(A, -b, 0.005, cg_max_iter*4, cgTFI_LogInterval*4)));
        res_norm_ratio = norm(g2c(P(:)).*dx(:))/norm(g2c(P(:)).*x_inner(:));
        x_inner = x_inner + dx;

        wres=m.*exp(1i*(real(ifftn(g2c(D).*fftn(g2c(P).*x_inner))))) - bi;

        cost_data_history(iter) = norm(wres(:),2);
        cost=abs(g2c(wG.*grad(P.*c2g(x_inner))));
        cost_reg_history(iter) = sum(cost(:));

        x(1:x_length) = reshape(x_inner,[x_length,1]);

        elapsed_time_sec = toc;
        minutes = floor(elapsed_time_sec / 60);
        seconds = floor(mod(elapsed_time_sec, 60));

        fprintf('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost:%8.4f\n elapsed_time: %d:%02d (min:sec).\n', ...
            iter, res_norm_ratio, cost_data_history(iter), cost_reg_history(iter), minutes, seconds);


    end
    if merit&&(iter<max_iter)
        wres = wres - mean(wres(mask(:)==1));
        a = wres(mask(:)==1);
        factor = std(abs(a))*6;
        wres = abs(wres)/factor;
        wres(wres<1) = 1;
        N_std(mask==1) = N_std(mask==1).*wres(mask==1).^2;
        tempn = double(N_std);
        m = dataterm_mask(data_weighting_mode, tempn, mask);
        b0 = m.*exp(1i*f);
    end
end

%%
x_map = P.*reshape(x(1:x_length),matrix_size)/(2*pi*delta_TE*CF)*1e6;
mask = logical(mask);
if flag_CSF
    x_map(mask) = x_map(mask) - mean(x_map(Mask_CSF));
end
if  ~flag_matrix_size
    x_map=x_map(1+size_change(1):matrix_size(1)-size_change(1), ...
        1+size_change(2):matrix_size(2)-size_change(2), ...
        1+size_change(3):matrix_size(3)-size_change(3));
    mask=mask(1+size_change(1):matrix_size(1)-size_change(1), ...
        1+size_change(2):matrix_size(2)-size_change(2), ...
        1+size_change(3):matrix_size(3)-size_change(3));
end
gTFI = x_map.*mask;

end