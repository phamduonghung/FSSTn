function fig4()
% fig4 : draws Figure 4 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.


clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

N = 1024;
t = (0:N-1)/N;

%% Tests signals
[a1,a2,if1,if2,s1,s2,~] = signal_test(t);
st1 = s1;

% Parameters
gamma = 0;
sigma = 0.05;
index = N/8+1:7*N/8;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Sparsity (no noise)
nc = 6;
es = 180;

s = st1;
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigma,ft,bt);
[~, ~, RM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(s,gamma,sigma);

plotsparse_new( RM(:,index),FSST2(:,index),FSST3(:,index),FSST4(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_1', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

%% 
s = s2;

[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigma,ft,bt);
[~, ~, RM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(s,gamma,sigma);

plotsparse_new( RM(:,index),FSST2(:,index),FSST3(:,index),FSST4(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_2', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

end
