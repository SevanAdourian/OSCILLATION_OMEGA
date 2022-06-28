clear; close all; clc

basename_re = "rho_ulm_re_lay";
basename_im = "rho_ulm_im_lay";

NR = 808; N_ICB = 166; N_CMB = 331;
del_rho_re_22 = zeros(NR,1);
del_rho_im_22 = zeros(NR,1);

prefactor = sqrt((5/4*pi)*(1/24));
% for i = 1:(NR-N_CMB)

for ii = 1:NR
    if (ii >N_CMB && ii < N_CMB+60)
        fname_re = strcat(strcat(basename_re,num2str(ii-(N_CMB+1),'%03.f')),".dat");
        fname_im = strcat(strcat(basename_im,num2str(ii-(N_CMB+1),'%03.f')),".dat");
    
        delta_rho_re = load(fname_re);
        delta_rho_im = load(fname_im);
        del_rho_re_22(ii) = delta_rho_re(3,3)*2;
        del_rho_im_22(ii) = delta_rho_im(3,3)*2;
    elseif (ii >= N_CMB+60)
        fname_re = strcat(strcat(basename_re,num2str(ii-(N_CMB+1),'%03.f')),".dat");
        fname_im = strcat(strcat(basename_im,num2str(ii-(N_CMB+1),'%03.f')),".dat");
    
        delta_rho_re = load(fname_re);
        delta_rho_im = load(fname_im);
        del_rho_re_22(ii) = delta_rho_re(3,3);
        del_rho_im_22(ii) = delta_rho_im(3,3);
    end
end
% for ii = 1:NR
%     if (ii >=N_CMB && ii <= N_CMB+30)
%     % if (ii >=N_CMB)
%         del_rho_re_22(ii) = (-0.016*0.1)/prefactor;
%         % del_rho_im_22(ii) = (-0.016*0.2)/prefactor;
%         del_rho_im_22(ii) = 0;
%     end
% end
% 
% dlmwrite('delta_rho_2_2_re_0_01_forte.dat',del_rho_re_22);
% dlmwrite('delta_rho_2_2_im_0_01_forte.dat',del_rho_im_22);


dlmwrite('delta_rho_2_2_s20rts_60_aug.dat',[del_rho_re_22, del_rho_im_22],'delimiter', '\t');
