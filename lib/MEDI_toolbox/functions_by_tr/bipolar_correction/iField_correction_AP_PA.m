filename=cell(0);
filename{1}= 'F:\MRIdata\TongRui\raw\meas_MID00202_FID35070_vibe_QSM_AP';
filename{2} = 'F:\MRIdata\TongRui\raw\meas_MID00203_FID35071_vibe_QSM_PA';
rawA_obj= mapVBVD(filename{1});
rawB_obj= mapVBVD(filename{2});
rawA = rawA_obj{1,2}.image(:,:,99,7,1);
rawB = rawB_obj{1,2}.image(:,:,99,7,1);
a = rawA(:,3);
b = rawB(:,3);
plot(abs(a),'b-');
hold on;
plot(abs(b),'r-');
a_fft = fftshift(fft(a));
b_fft = fftshift(fft(b));
plot(abs(a_fft));
hold off;
plot(abs(b_fft));
c = angle(a./b);
% 
% [m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
% c(ind==1)=c(ind==1)-2*pi;
% c(ind==3)=c(ind==3)+2*pi;
% 
% for n=1:min(2,nechos-1)
%     cd=((Y(:,n+1)-Y(:,n)))-c;
%     Y(cd<-pi,(n+1):end)=Y(cd<-pi,n+1:end)+2*pi;
%     Y(cd>pi,(n+1):end)=Y(cd>pi,n+1:end)-2*pi;
% end
figure(1)
for i=1:20
a = rawA_obj{1,2}.image(:,i,:,:,1);
a_fft = fftshift(fftn(a));
a_fft = squeeze(a_fft);
subplot(2,10,i)
imshow(abs(a_fft(:,:,5)));
end

figure(2)
for i=1:20
b = rawB_obj{1,2}.image(:,i,:,:,1);
b_fft = fftshift(fftn(b));
b_fft = squeeze(b_fft);
subplot(2,10,i)
imshow(abs(b_fft(:,:,5)));
end


