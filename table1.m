function table1()

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
N = 1024;
gamma = 0;
sigma = 0.05;
index = N/8+1:7*N/8;
d = 0;
clwin = 10;
lambda = 0;
nmodes = 2;

t = (0:N-1)/N;
% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
 [a1,a2,if1,if2,s1,s2,s] = signal_test(t);

%% TF presentations of signal
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigma,ft,bt);

%% Ridge extraction FSST2

[Cs, Es] = exridge_mult(FSST2, nmodes, lambda, clwin);
imf2 = sigma*recmodes(FSST2,Cs,d);

disp('SNRout by FSST2')
disp(snr(real(s1(index)),real(imf2(2,index)- s1(index))));
disp(snr(real(s2(index)),real(imf2(1,index)- s2(index))));
disp(snr(real(s(index)),real(imf2(2,index)+imf2(1,index)- s(index))));


%% Ridge extraction FSST3

[Cs, Es] = exridge_mult(FSST3, nmodes, lambda, clwin);
imf3 = sigma*recmodes(FSST3,Cs,d);

disp('SNRout by FSST3')
disp(snr(real(s1(index)),real(imf3(2,index)- s1(index))));
disp(snr(real(s2(index)),real(imf3(1,index)- s2(index))));
disp(snr(real(s(index)),real(imf3(2,index)+imf3(1,index)- s(index))));

%% Ridge extraction FSST4

[Cs, Es] = exridge_mult(FSST4, nmodes, lambda, clwin);
imf4 = sigma*recmodes(FSST4,Cs,d);

disp('SNRout by FSST4')
disp(snr(real(s1(index)),real(imf4(2,index)- s1(index))));
disp(snr(real(s2(index)),real(imf4(1,index)- s2(index))));
disp(snr(real(s(index)),real(imf4(2,index)+imf4(1,index)- s(index))));

%% Mode reconstruction accuracy
ac = zeros(3,3);
%FSST2
ac(1,1) = snr(real(s1(index)),real(imf2(2,index)- s1(index)));
ac(2,1) = snr(real(s2(index)),real(imf2(1,index)- s2(index)));
ac(3,1) = snr(real(s(index)),real(imf2(2,index)+imf2(1,index)- s(index)));
%FSST3
ac(1,2) = snr(real(s1(index)),real(imf3(2,index)- s1(index)));
ac(2,2) =  snr(real(s2(index)),real(imf3(1,index)- s2(index)));
ac(3,2) = snr(real(s(index)),real(imf3(2,index)+imf3(1,index)- s(index)));
%FSST4
ac(1,3) = snr(real(s1(index)),real(imf4(2,index)- s1(index)));
ac(2,3) =  snr(real(s2(index)),real(imf4(1,index)- s2(index)));
ac(3,3) = snr(real(s(index)),real(imf4(2,index)+imf4(1,index)- s(index)));

digits(3); %this changes the output precision
latex(sym(ac,'d'))
end
