function [water fat freq R2s ErrorMap iter model]= fit_IDEAL_R2_bip(s0, t,CF, f0, R2s, max_iter,ErrorMap)
matrix_size = size(s0);
numvox = prod(matrix_size(1:end-1));
numte = matrix_size(end);

if nargin<7
    ErrorMap = zeros([1 numvox]);
end

if nargin<6
    max_iter =5;
end
if nargin <5
    R2s = zeros([1 numvox]);
end
if nargin<4 | isempty(f0)
    f0 = zeros([1 numvox]);
end
if nargin<3
    f_fat = -420;
end
    
if numel(f0) == 1
    f0 = f0*ones([1 numvox]);
end
if numel(R2s) == 1
    R2s = R2s*ones([1 numvox]);
end
if numel( ErrorMap) == 1
     ErrorMap =  ErrorMap*ones([1 numvox]);
end

s0 = permute(s0, [length(matrix_size) 1:length(matrix_size)-1]);
s0 = reshape(s0,[numte numvox]);
R2s = reshape(R2s,[1 numvox]);
f0 = reshape(f0,[1 numvox]);
ErrorMap(isnan(ErrorMap)) = 0;
ErrorMap = reshape(ErrorMap,[1 numvox]);
t = reshape(t,[numte 1]);
t = repmat(t,[1 numvox]);


EchoNum = (1:numte)';
EchoNum = repmat(EchoNum,[1, numvox]);

% O = ones([numte numvox]).*exp(repmat(-R2s,[numte 1]).*t);
% C = real(exp(-1i*2*pi*f_fat*t)).*exp(repmat(-R2s, [numte 1]).*t);

O = ones([numte numvox]);
%µ¥·å
% f_fat = -3.5*1e-6*CF;
% C = exp(-1i*2*pi*f_fat*t);

%6·å
% fat_chemshift = [0.90, 1.30, 2.10, 2.76, 4.31, 5.30]-4.7;
% fat_frequency=CF*1e-6*fat_chemshift;
% fat_relAmps = [0.087, 0.693, 0.128, 0.004, 0.039, 0.048];
% C=0;
% for k = 1:6
%    C = C + fat_relAmps(k).* exp(-1i*2*pi*fat_frequency(k).*t);
% end;

%7·å
%  fat_chemshift = [0.90, 1.30, 1.59, 2.03, 2.25, 2.77, 5.31] -4.7;
%  fat_frequency=CF*1e-6*fat_chemshift ;
%  fat_relAmps = [0.083, 0.627, 0.072, 0.096, 0.066, 0.015, 0.042];
%  C=0;
% for k = 1:7
%    C = C + fat_relAmps(k).* exp(-1i*2*pi*fat_frequency(k).*t);
% end;

%9·å
fat_chemshift = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29] - 4.7;
  fat_frequency=CF*1e-6*fat_chemshift ;
fat_relAmps = [0.088, 0.642, 0.058, 0.062, 0.058, 0.006, 0.039, 0.01, 0.037];
C=0;
for k = 1:9
   C = C + fat_relAmps(k).* exp(-1i*2*pi*fat_frequency(k).*t);
end;


y = zeros([4 numvox]);
dy = zeros([4 numvox]);
dy(1,:) = 1e4;
iter = 0;

y(1,:) = f0+1i*R2s/(2*pi);
y(2,:) = ErrorMap;

%E =1;
E = exp(1i *(-1).^EchoNum .*repmat(y(2,:),[numte,1]));         %phaseErrorMap
P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
y(3:4,:) = invA(P.*E.*O, P.*E.*C, s0); 

% PP = diag(exp(-1i*2*pi*(f0+1i*R2s/(2*pi))*tt))
% 
% AA = [ones(6,1),exp(-1i*2*pi*f_fat*tt)]
% M = PP*AA;
% N = inv(M'*M)*M'*s0;

update = dy(1,:);

while(iter<max_iter)
%while (iter<max_iter)&&(sqrt(sum(real(update).^2,2)/numvox)>0.1)
    sn = P.*E.*O.*repmat(y(3,:),[numte 1]) + P.*E.*C.*repmat(y(4,:),[numte 1]);
    sr = s0 - sn;

%     gr = -2*pi*t.*(-repmat(y(2,:),[numte 1]).*z - repmat(y(3,:),[numte 1]).*O - repmat(y(4,:),[numte 1]).*d - repmat(y(5,:),[numte 1]).*C);
%     gi = -2*pi*t.*(repmat(y(2,:),[numte 1]).*O - repmat(y(3,:),[numte 1]).*z + repmat(y(4,:),[numte 1]).*C - repmat(y(5,:),[numte 1]).*d);

    Bcol01 = -1i*2*pi*t.*sn;
    Bcol02 = 1i*(-1).^EchoNum.*sn;
    dy = invC(Bcol01,Bcol02, P.*E.*O, P.*E.*C, sr);     
    y = y+dy;
    iter = iter+1;

    temp = y(1,:);
    temp(abs(imag(temp))*2*pi>50)=real(temp(abs(imag(temp))*2*pi>50));
%     y(1,:) = temp;  hann_low
    
   % y(2,:)=1i*abs(exp(1i.*y(2,:)))+angle(exp(1i.*y(2,:)));
    y(2,:)=angle(exp(1i.*y(2,:)));
    E = exp(1i *(-1).^EchoNum .*repmat(y(2,:),[numte,1]));         %phaseErrorMap
    P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
    y(3:4,:) = invA(P.*E.*O, P.*E.*C, s0);

    update = dy(1,:);
    update(isnan(update)) = 0;
    update(isinf(update)) = 0;
end

freq = reshape(real(y(1,:)),matrix_size(1:end-1));
R2s = -reshape(imag(y(1,:))*2*pi,matrix_size(1:end-1));
water = reshape(y(3,:),matrix_size(1:end-1));
fat = reshape(y(4,:),matrix_size(1:end-1));
ErrorMap = reshape(y(2,:),matrix_size(1:end-1));
model = P.*E.*(O.*repmat(y(3,:),[numte 1]) + C.*repmat(y(4,:),[numte 1]));
model = reshape(model,matrix_size);
% s_model = [s0;model];
% figure; plot(s0,'ro'); hold on; plot(model,'bx');hold off;
% axis([-max(abs(real(s_model))) max(abs(real(s_model))) -max(abs(imag(s_model))) max(abs(imag(s_model)))]*1.2)

freq(isinf(freq)) = 0;
freq(isnan(freq)) = 0;
freq(abs(freq)>10e4) = 0;
R2s(isinf(R2s)) = 0;
R2s(isnan(R2s)) = 0;
water(isinf(water)) = 0;
water(isnan(water)) = 0;
water(abs(water)>10e5) = 0;
fat(isinf(fat)) = 0;
fat(isnan(fat)) = 0;
fat(abs(fat)>10e5) = 0;
ErrorMap(isinf(ErrorMap))=0;
ErrorMap(isnan(ErrorMap)) =0;
model(isinf(model)) = 0;
model(isnan(model)) = 0;

% freq(isinf(freq)) = 0;
% freq(isnan(freq)) = 0;
% 
% R2s(isinf(R2s)) = 0;
% R2s(isnan(R2s)) = 0;
% water(isinf(water)) = 0;
% water(isnan(water)) = 0;
% 
% fat(isinf(fat)) = 0;
% fat(isnan(fat)) = 0;
% 
% model(isinf(model)) = 0;
% model(isnan(model)) = 0;


function x=invA(col1, col2, y)
% assemble A^H*A
a11 = sum(conj(col1).*col1, 1);
a12 = sum(conj(col1).*col2, 1);
a22 = sum(conj(col2).*col2, 1);

% inversion of A^H*A
d = (a11.*a22 - a12.*conj(a12));
ia11 = a22./d;
ia12 = -a12./d;
ia22 = a11./d;

% y project onto A^H
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
% calculate x
x(1,:) = sum(ia11.*py1 + ia12.*py2, 1);
x(2,:) = sum(conj(ia12).*py1 + ia22.*py2, 1);

function x=invB(col1, col2, col3, y)
% assemble B^H*B
b11 = sum(conj(col1).*col1,1);
b12 = sum(conj(col1).*col2,1);
b13 = sum(conj(col1).*col3,1);
b22 = sum(conj(col2).*col2,1);
b23 = sum(conj(col2).*col3,1);
b33 = sum(conj(col3).*col3,1);

% inversion of B'*B
d = (b13.*conj(b12).*conj(b23) + b11.*b22.*b33 + b12.*b23.*conj(b13) - b13.*b22.*conj(b13) - b11.*b23.*conj(b23) - b12.*b33.*conj(b12));
ib11 = (b22.*b33 - b23.*conj(b23))./d;
ib12 = -(b12.*b33 - b13.*conj(b23))./d;
ib13 = (b12.*b23 - b13.*b22)./d;
ib22 = (b11.*b33 - b13.*conj(b13))./d;
ib23 = -(b11.*b23 - b13.*conj(b12))./d;
ib33 =  (b11.*b22 - b12.*conj(b12))./d;

% y project onto B'
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
py3 = sum(conj(col3).*y,1);
% calculate x
x(1,:) = sum(ib11.*py1 + ib12.*py2 + ib13.*py3, 1);
x(2,:) = sum(conj(ib12).*py1 + ib22.*py2 + ib23.*py3, 1);
x(3,:) = sum(conj(ib13).*py1 + conj(ib23).*py2 + ib33.*py3, 1);

function x=invC(col1, col2, col3, col4, y)
% assemble A'*A
a11 = sum(col1.*col1,1);
a12 = sum(col1.*col2,1);
a13 = sum(col1.*col3,1);
a14 = sum(col1.*col4,1);
a22 = sum(col2.*col2,1);
a23 = sum(col2.*col3,1);
a24 = sum(col2.*col4,1);
a33 = sum(col3.*col3,1);
a34 = sum(col3.*col4,1);
a44 = sum(col4.*col4,1);

% inversion of A'*A
d = (a33.*a44.*a12.^2 - a12.^2.*a34.^2 - 2.*a44.*a12.*a13.*a23 + 2.*a12.*a13.*a24.*a34 + 2.*a12.*a14.*a23.*a34 - 2.*a33.*a12.*a14.*a24 - a13.^2.*a24.^2 + a22.*a44.*a13.^2 + 2.*a13.*a14.*a23.*a24 - 2.*a22.*a13.*a14.*a34 - a14.^2.*a23.^2 + a22.*a33.*a14.^2 + a11.*a44.*a23.^2 - 2.*a11.*a23.*a24.*a34 + a11.*a33.*a24.^2 + a11.*a22.*a34.^2 - a11.*a22.*a33.*a44);
ia11 = (a44.*a23.^2 - 2.*a23.*a24.*a34 + a33.*a24.^2 + a22.*a34.^2 - a22.*a33.*a44)./d;
ia12 = -(a12.*a34.^2 - a13.*a24.*a34 - a14.*a23.*a34 + a14.*a24.*a33 + a13.*a23.*a44 - a12.*a33.*a44)./d;
ia13 = -(a13.*a24.^2 - a14.*a23.*a24 - a12.*a24.*a34 + a14.*a22.*a34 + a12.*a23.*a44 - a13.*a22.*a44)./d;
ia14 = -(a14.*a23.^2 - a13.*a23.*a24 - a12.*a23.*a34 + a12.*a24.*a33 + a13.*a22.*a34 - a14.*a22.*a33)./d;
ia22 = (a44.*a13.^2 - 2.*a13.*a14.*a34 + a33.*a14.^2 + a11.*a34.^2 - a11.*a33.*a44)./d;
ia23 = -(a14.^2.*a23 - a13.*a14.*a24 - a12.*a14.*a34 + a11.*a24.*a34 + a12.*a13.*a44 - a11.*a23.*a44)./d;
ia24 = -(a13.^2.*a24 - a13.*a14.*a23 - a12.*a13.*a34 + a12.*a14.*a33 + a11.*a23.*a34 - a11.*a24.*a33)./d;
ia33 = (a44.*a12.^2 - 2.*a12.*a14.*a24 + a22.*a14.^2 + a11.*a24.^2 - a11.*a22.*a44)./d;
ia34 = -(a12.^2.*a34 - a12.*a13.*a24 - a12.*a14.*a23 + a13.*a14.*a22 + a11.*a23.*a24 - a11.*a22.*a34)./d;
ia44 = (a33.*a12.^2 - 2.*a12.*a13.*a23 + a22.*a13.^2 + a11.*a23.^2 - a11.*a22.*a33)./d;
% y project onto A'
py1 = sum(col1.*y,1);
py2 = sum(col2.*y,1);
py3 = sum(col3.*y,1);
py4 = sum(col4.*y,1);
% calculate x
x(1,:) = sum(ia11.*py1 + ia12.*py2 + ia13.*py3 + ia14.*py4,1);
x(2,:) = sum(ia12.*py1 + ia22.*py2 + ia23.*py3 + ia24.*py4,1);
x(3,:) = sum(ia13.*py1 + ia23.*py2 + ia33.*py3 + ia34.*py4,1);
x(4,:) = sum(ia14.*py1 + ia24.*py2 + ia34.*py3 + ia44.*py4,1);
