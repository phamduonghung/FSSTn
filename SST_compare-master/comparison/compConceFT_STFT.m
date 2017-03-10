% This code draws the ConceFT transform based on the synchrosqueezed short-time Fourier transform

L = 1 ;
load dat.mat;
N = length(x) ;
Hz = N ;
highFq = 130;
WinLen = 401;%301 ;
        %% Setup parameters
        %% alpha is the grid size in the freq axis
alpha = 0.025/Hz ;
xm = x ;
time = [1/Hz:1/Hz:L]' ;
MT = 10; %MT = 1: ordinary SST; MT > 1: ConceFT % number of random projections
hop = 1;
dim = 2; %1; 5% maximum order % number of orthogonal mother wavelet
supp = 6; %half time support (>= 6 recommended) 
Smooth = 1; % 0 or 1
Hemi = 1; % 0 or 1

%%%% call sqSTFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(xm,  0, 0.5*highFq/Hz*2, alpha, hop, WinLen, dim, supp, MT, Smooth, Hemi) ;
pic = figure;
imageSQ(time, tfrsqtic*Hz, abs(ConceFT)) ; axis square; 
xlabel('Time') ; ylabel('Freq') ;
title('MT=10,dim=2');
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsSTFT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);


MT = 30; %MT = 1: ordinary SST; MT > 1: ConceFT % number of random projections
hop = 1;
dim = 2; %1; 5% maximum order % number of orthogonal mother wavelet
supp = 6; %half time support (>= 6 recommended) 
Smooth = 1; % 0 or 1
Hemi = 1; % 0 or 1

%%%% call sqSTFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(xm,  0, 0.5*highFq/Hz*2, alpha, hop, WinLen, dim, supp, MT, Smooth, Hemi) ;
pic = figure;
imageSQ(time, tfrsqtic*Hz, abs(ConceFT)) ; axis square; 
xlabel('Time') ; ylabel('Freq') ;
title('MT=30,dim=2');
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsSTFT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);


