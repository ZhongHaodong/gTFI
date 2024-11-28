function rho = Search_rho_plus( ibeta,middle_slice, D, F )
% 互相关函数
        
%         a_kspace = D{middle_slice-5:middle_slice+5,1};
%         b_kspace = D{middle_slice-5:middle_slice+5,2};
%         c_kspace = D{middle_slice-5:middle_slice+5,3};
middle_slice = middle_slice -10;
    rho = 0;
    alpha =0*pi;      
    beta = ibeta*2*pi / size(F,1);
    rows = int16(size(F,1)/2) -  int16(size(F,1)/4)  :  int16(size(F,1)/2) + int16(size(F,1)/4) ;
    phase = @(f,alpha,beta)              exp(1i* (alpha + beta * f));          %%附加相位函数
    for i = 0:3
        a_kspace = D{middle_slice +i,1};
        b_kspace = D{middle_slice +i,2};
        c_kspace = D{middle_slice +i,3};               
        E2= ifft(fft(b_kspace).*phase(F,alpha,beta));
        % tmp1 = abs(corr(a_kspace(90:140,60:100),E2(90:140,60:100)));
        % tmp2 = abs(corr(b_kspace(90:140,60:100),E2(90:140,60:100)));         
         tmp1 = abs(corr(a_kspace(rows,:),E2(rows,:)));
         tmp2 = abs(corr(c_kspace(rows,:),E2(rows,:)));
         rho = -sum(diag(tmp1).*diag(tmp2)) + rho;
    end

end


