% This code draws the synchrosqueezed short-time Fourier transform in Figure 1 of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

L = 1 ;
load dat.mat;
N = length(x) ;
Hz = N ;
highFq = 130;
opts.WinLen = 401;%301 ;
        %% Setup parameters
        %% alpha is the grid size in the freq axis
alpha = 0.025/Hz ;
xm = x ;
time = [1/Hz:1/Hz:L]' ;


%%%% call sqSTFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h,Dh,t] = hermf(opts.WinLen, 1, 6) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(xm, 0, 0.5*highFq/Hz*2, alpha, 1, h', Dh') ; 
figure;
imageSQ(time, tfrsqtic*Hz, abs(tfrsq)) ; axis square; 
xlabel('Time') ; ylabel('Freq') ;
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsSTFT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

