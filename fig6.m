function fig6()
% fig6 : draws Figure 6 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);

N = 1024;

chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
gamma = 10^(-2);
sigma = 0.05;
index = N/8+1:7*N/8;
d = 0:1:8;
clwin = 10;
lambda = 0;
nmodes = 2;

t = (0:N-1)/N;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 2
[a1,a2,if1,if2,s1,s2,st1] = signal_test(t);
s=st1;

%% TF presentations of signal
[~,~,FSST2,FSST3,FSST4,~,~,~,~,~,~,~,~] = sstn(s,gamma,sigma,ft,bt);

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

SNR_FSST2 = zeros(length(d),1);
SNR_FSST3 = zeros(length(d),1);
SNR_FSST4 = zeros(length(d),1);
[Cs2, ~] = exridge_mult(FSST2, nmodes, lambda, clwin);
[Cs3, ~] = exridge_mult(FSST3, nmodes, lambda, clwin);
[Cs4, ~] = exridge_mult(FSST4, nmodes, lambda, clwin);
close all;

for k=1:length(d)
    imf2 = sigma*recmodes(FSST2,Cs2,d(k));
    SNR_FSST2(k) = snr(real(s(index)),real(imf2(1,index)+imf2(2,index)- s(index)));
     
    imf3 = sigma*recmodes(FSST3,Cs3,d(k));
    SNR_FSST3(k) = snr(real(s(index)),real(imf3(1,index)+imf3(2,index)- s(index)));

    imf4 = sigma*recmodes(FSST4,Cs4,d(k));
    SNR_FSST4(k) = snr(real(s(index)),real(imf4(1,index)+imf4(2,index)- s(index)));
end

FigHandle = figure; 
%set(FigHandle, 'Position', [400,400, 800, 600]);

%% SNR_out

plot(d,SNR_FSST2,'bs-',d,SNR_FSST3,'gd-',d,SNR_FSST4,'rd-','LineWidth',1);

legend('FSST2','FSST3','FSST4','Location','northwest'); 
ylabel('output SNR (dB)');
xlabel('d');

%%%%%%%%%%%%%%%%%%%%%% print SNR out signal
export_fig(FigHandle, ... % figure handle
    sprintf('%s/SNRout_d', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
end
