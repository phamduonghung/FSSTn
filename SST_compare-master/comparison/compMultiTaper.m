% This code draws the multitaper (reassigned) spectrograms in Figure 8 of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.
%
% This code adapts the demo code by
% P. Flandrin & J. Xiao, 2005
%
% computes and displays Hermite multitaper (reassigned) spectrograms
% in the case of a synthetic sinusoidal FM signal embedded in white
% Gaussian noise
%

%set up data
load dat.mat;
Nx = length(f2);
t = 1:Nx;
signal = f2.';

% analysis parameters

Nfft = Nx*5 ; % number of FFT bins
M = 10 ; % number of Hermite tapers
Nh = 505; % length of Hermite functions

% computation
tic;
[S,RS,hat,tt] = tfrrsp_hm(x,t,Nfft,Nh,M,6) ;

SS = [] ; RSS = [] ;


for opt = 1:2 % loop on means (cf. mean_hm.m)
    for k = M % loop on tapers
        [opt k]
        [Sm,RSm] = mean_hm(S(:,:,1:k),RS(:,:,1:k),hat(:,:,1:k),opt) ;
        SS(:,:,opt,k) = Sm ;
        RSS(:,:,opt,k) = RSm ;
    end
end
toc
% display

col = 1 ;
for opt = 1:2 % loop on means
    for k = M % loop on # tapers
        figure;
        imagesc([0,1],[0,130],RSS(1:Nfft/8,:,opt,k));axis square; axis xy;
        xlabel('Time')
        ylabel('Freq')
%         set(gca, 'FontSize', 18);
%         b=get(gca);
%         set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
%         str = ['results/nsRS' num2str(opt)];
%         print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    end
end

