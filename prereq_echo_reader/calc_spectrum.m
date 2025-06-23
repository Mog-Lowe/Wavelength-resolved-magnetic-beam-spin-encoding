function [fx,fy]=calc_spectrum(curr,beff,signal,molec,T)
 load_chess_parameters;
 if strcmpi(molec,'H2')
     mass=2*SE_amu* 1.00794;
%      display('PARAMETES OF H2')
 elseif strcmpi(molec,'He3')
     mass=SE_3hemass;
%      display('PARAMETES OF He3')
 elseif strcmpi(molec,'D2')
     mass=2*SE_amu*2.01410178;
%      display('PARAMETES OF D2')
 elseif strcmpi(molec,'He4')
     mass=SE_amu*4.002602
%      display('PARAMETES OF He4')
 end
curr=curr-min(curr);
B=curr*beff;     
dB=B(2)-B(1);
fB=max(B)-min(B); 
% T=value;

if mass==SE_3hemass
    gamma=32.43*10^6;
    spec=abs(fftshift(fft(signal-mean(signal))));
    dlambda=SE_h/(gamma*mass*fB);
    [zrv,zrp]=min(spec);
    spec_pos=spec(zrp:end);
    lambda_pos=(0:length(spec_pos)-1)*dlambda;
    fy=spec_pos;
    fx=lambda_pos;
else
    Ei_SI=5/2*SE_kB*T;
    ki_SI = (2*mass*Ei_SI).^0.5/SE_hbar;
    lambda=2*pi/(ki_SI);
    spec=abs(fftshift(fft(signal-mean(signal))));
    dgamma=SE_h/(lambda*mass*fB);
    [zrv,zrp]=min(spec);
    spec_pos=spec(zrp:end);
    gamma_pos=(0:length(spec_pos)-1)*dgamma;
    fy=spec_pos;
    fx=gamma_pos;
end
end