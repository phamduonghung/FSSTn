% fig5 : draws Figure 5 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.
%
% It compares the ideal TF representation with the one given by different
% time-frequency representations, using the earth-mover distance.
% 
% This script needs some functions, including the fast EMD and the
% synchrosqueezed wavelet packet transform, that can be downloaded from 
% https://github.com/HaizhaoYang/SST_compare/blob/master/comparison/

function fig5()

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

%% set up data
N  = 1024;
t  = (0:N-1)/N;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;
 
 %% Test signal 
[a1,a2,if1,if2,s1,s2,sss] = signal_test(t);
s = real(s1);
InstFreq = if1;
i=1;
% simulation parameters
numTest = 1; % Nb realizations 
NMvec = -5:5:30;
index = N/8+1:7*N/8;

% store results
resFSST2 = zeros(length(NMvec),1);
resFSST3 = zeros(length(NMvec),1);
resFSST4 = zeros(length(NMvec),1);
resRM = zeros(length(NMvec),1);
resSSWPT1 = zeros(length(NMvec),1);
resSSWPT20 = zeros(length(NMvec),1);

%% parameters for SSWPT
res = 1;
NG = N;
is_unif = 1;
typeNUFFT = 1;
is_cos = 1;
epsl = 1e-2;
xo = s;
rad = 1; % 1.3
fff = s.';
red = 1; % 10
t_sc = 1/2 + 2/8;
lowfreq = 0; highfreq = N/2;
is_real = 1;

for cntNM = 1:length(NMvec) 
    % loop on noise
    
     sg = sqrt( var(s)*10^(-NMvec(cntNM)/10));
    for cntt = 1:numTest
        % different realizations
        
        % set noise
        b = sg*randn(size(s));
        sb = s+b;
        SNR = round(10*log10((var(s(:)))/var(b(:))))
        
        %% Computes FSST and FSST2
        % parameters
        gamma = 0;
        sigma = 0.05;

        % FSST2, FSST3, FSST4
        [~,~,FSST2,FSST3,FSST4,~,~,~,~,~,~,~,~] = sstn(sb,gamma,sigma,ft,bt);
        [~, ~, RM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(sb,gamma,sigma);

        % SSWPT
        SSWPT1 = ss_wp1_fwd(real(sb),is_real,is_unif,typeNUFFT,t,N,highfreq,lowfreq,rad,is_cos,t_sc,1,epsl,res,0);
        SSWPT20 = ss_wp1_fwd(real(sb),is_real,is_unif,typeNUFFT,t,N,highfreq,lowfreq,rad,is_cos,t_sc,20,epsl,res,0);
                
        % computes EMD
        resFSST2(cntNM) = resFSST2(cntNM) + EMDMat(abs(FSST2(:,index)),InstFreq(index),N/2)/numTest;
        resFSST3(cntNM) = resFSST3(cntNM) + EMDMat(abs(FSST3(:,index)),InstFreq(index),N/2)/numTest;
        resFSST4(cntNM) = resFSST4(cntNM) + EMDMat(abs(FSST4(:,index)),InstFreq(index),N/2)/numTest;
        resSSWPT1(cntNM) = resSSWPT1(cntNM) + EMDMat(abs(SSWPT1(:,index)),InstFreq(index),N/2)/numTest;
        resSSWPT20(cntNM) = resSSWPT20(cntNM) + EMDMat(abs(SSWPT20(:,index)),InstFreq(index),N/2)/numTest;
        resRM(cntNM) = resRM(cntNM) + EMDMat(abs(RM(:,index)),InstFreq(index),N/2)/numTest;
        
        % inspects results
        if 0
            figure;
            subplot(221);imagesc(log(0.1+abs(FSST2(:,index))));title('FSST2');
            subplot(222);imagesc(log(0.1+abs(FSST3(:,index))));title('FSST3');
            subplot(223);imagesc(log(0.1+abs(FSST4(:,index))));title('FSST4');
            subplot(224);imagesc(log(0.1+abs(resRM(:,index))));title('RM');
        end

    end
end

set(findall(0,'type','axes'),'xtick',[],'ytick',[]);

FigHandle7 = figure();
%set(FigHandle7, 'Position', [400,400, 800, 450]);

plot(NMvec,resRM,'m-d',NMvec,resFSST2,'k-*',NMvec,resFSST3,'b-o',NMvec,resFSST4,'r-+');
legend('RM','FSST2','FSST3','FSST4','Location','southwest'); xlim([0, 30]);
xlabel('input SNR');ylabel('EMD'); %title('(b)')

%%%%%%%%%%%%%%%%%%%%%%print EDM
export_fig(FigHandle7, ... % figure handle
    sprintf('%s/EMD_%d', chemin0,i),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
end
