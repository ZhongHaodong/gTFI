function gTFI = gTFI_matrixForm(params_gTFI_MatrixForm)
mesh = params_gTFI_MatrixForm.mesh;

delta_TE = params_gTFI_MatrixForm.delta_TE;
CF = params_gTFI_MatrixForm.CF;
voxel_size = params_gTFI_MatrixForm.voxel_size;
mask = logical(params_gTFI_MatrixForm.mask);
iMag = params_gTFI_MatrixForm.iMag;
B0_dir = params_gTFI_MatrixForm.B0_dir;
L = params_gTFI_MatrixForm.L;
%% INVERSE
%%%% initialization %%%%%%%%
%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
RDF = params_gTFI_MatrixForm.f;
matrix_size = size(RDF);
lambda = params_gTFI_MatrixForm.lambda;

cg_max_iter = getOrDefault(params_gTFI_MatrixForm, 'cg_max_iter', 500);
cg_tol = getOrDefault(params_gTFI_MatrixForm, 'cg_tol', 0.01);
max_iter = getOrDefault(params_gTFI_MatrixForm, 'max_iter', 6);
tol_norm_ratio = getOrDefault(params_gTFI_MatrixForm, 'tol_norm_ratio', 0.05);

data_weighting_mode = 1;
grad = @fgrad;
div = @bdiv;

N_std = params_gTFI_MatrixForm.N_std;
tempn = double(N_std);
D=dipole_kernel(matrix_size, voxel_size, B0_dir);
m = dataterm_mask(data_weighting_mode, tempn, mask);

iter=0;
x_length = matrix_size(1)*matrix_size(2)*matrix_size(3);

f_length = size(mesh.faces,1);
t_length = x_length+ f_length;
x = zeros([t_length ,1]); %real(ifftn(conj(D).*fftn((abs(m).^2).*RDF)));
zero = zeros([f_length,1]);

b0 = m.*exp(1i*RDF);

res_norm_ratio_x = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);

e=0.000001; %a very small number to avoid /0
Dconv = @(dx) real(ifftn(D.*fftn(reshape(dx(1:x_length),matrix_size))))+u_map(L*dx(x_length+1:end),mask);
DconvH = @(dx) [reshape(real(ifftn(D.*fftn(dx))),[x_length,1]);L'*dx(mask)];

%初始化边界条件
A= @(x) L'*u_flat(mask.*u_map(L*x,mask),mask);
x(x_length+1:end) = real(pcg(A,L'*u_flat(mask.*(RDF-u_map(L*x(x_length+1:end),mask)),mask),0.01,50,[],[]));

wG = gradient_mask(1, iMag, mask, grad, voxel_size);
while (res_norm_ratio_x>tol_norm_ratio)&&(iter<max_iter)
    tic
    iter=iter+1;
    Vr = 1./sqrt(abs(wG.*grad(real(reshape(x(1:x_length),matrix_size)),voxel_size)).^2+e);
    w = m.*exp(1i*(ifftn(D.*fftn(reshape(x(1:x_length),matrix_size)))+u_map(L*x(x_length+1:end),mask)));
    reg = @(dx) [reshape(div(wG.*(Vr.*(wG.*grad(real(reshape(dx(1:x_length),matrix_size)),voxel_size))),voxel_size),[x_length,1]);zero];

    fidelity = @(dx)DconvH(conj(w).*w.*Dconv(dx) );

    A =  @(dx) reg(dx) + 2*lambda*fidelity(dx);
    b = reg(x) + 2*lambda*DconvH( real(conj(w).*conj(1i).*(w-b0)) );

    dx = real(cgsolve(A, -b, cg_tol, cg_max_iter, 100));
    x = x + dx;
    
    dx2 = [reshape(reshape(dx(1:x_length),matrix_size),[x_length,1]);dx(x_length+1:end)];
    x2 = [reshape(reshape(x(1:x_length),matrix_size),[x_length,1]);x(x_length+1:end)];
    res_norm_ratio = norm(dx2(:))/norm(x2(:));
    res_norm_ratio_x = norm(dx2(1:x_length))/norm(x2(1:x_length));
    res_norm_ratio_b = norm(dx2(x_length:end))/norm(x2(x_length:end));


    wres=m.*exp(1i*(real(ifftn(D.*fftn(reshape(x(1:x_length),matrix_size))))+u_map(L*x(x_length+1:end),mask))) - b0;

    cost_data_history(iter) = norm(wres(:),2);
    cost=abs(grad(x));
    cost_reg_history(iter) = sum(cost(:));

    if 1
        wres = wres - mean(wres(mask(:)==1));
        a = wres(mask(:)==1);
        factor = std(abs(a))*6;
        wres = abs(wres)/factor;
        wres(wres<1) = 1;
        N_std(mask==1) = N_std(mask==1).*wres(mask==1).^2;
        tempn = double(N_std);
        m = dataterm_mask(data_weighting_mode, tempn, mask);
        b0 = m.*exp(1i*RDF);
    end

    fprintf('iter: %d;\n res_norm_ratio:%8.4f; res_norm_ratio_x:%8.4f; res_norm_ratio_b:%8.4f;\n cost_L2:%8.4f; cost:%8.4f.\n'...
            ,iter, res_norm_ratio, res_norm_ratio_x, res_norm_ratio_b,cost_data_history(iter), cost_reg_history(iter));
end
%convert x to ppm
gTFI = reshape(x(1:x_length),matrix_size)/(2*pi*delta_TE*CF)*1e6.*mask;

end
