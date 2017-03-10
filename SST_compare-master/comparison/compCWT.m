% This code draws the synchrosqueezed wavelet transform in Figure 1 of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

opts = struct();
opts.motherwavelet = 'Cinfc' ;
opts.CENTER = 1 ;
opts.FWHM = 0.3 ;

lowfreq = 0 ;
highfreq = 130 ;

%% 
Hz = 2^10 ;
L = 1 ;
time = [1/Hz:1/Hz:L]' ;
N = length(time) ;
alpha = 0.2 ; 

scrsz = get(0,'ScreenSize');

load dat.mat;
xm = x;

[~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
figure;
imageSQ(time, tfrsqtic, abs(tfrsq)); axis square;
xlabel('Time') ; ylabel('Freq') ; 
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsCWT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
