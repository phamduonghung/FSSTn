% This code draws the ConceFT transform based on the synchrosqueezed wavelet transform

opts = struct();
opts.motherwavelet = 'Cinfc' ;
opts.CENTER = 1 ;
opts.FWHM = 0.3 ;
opts.dim = 2;%1; 5% maximum order

lowfreq = 0 ;
highfreq = 130 ;

%% 
Hz = 2^10 ;
L = 1 ;
time = [1/Hz:1/Hz:L]' ;
N = length(time) ;
alpha = 0.2 ;
MT = 10;
Smooth = 1;
Hemi = 1;

scrsz = get(0,'ScreenSize');

load dat.mat;
xm = x;

[tfrsq, ConceFT, tfrsqtic] = ConceFT_CWT(time, xm, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi) ;
%[~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
figure;
imageSQ(time, tfrsqtic, abs(ConceFT)); axis square;
xlabel('Time') ; ylabel('Freq') ; 
title('MT=10');
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsCWT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

MT = 30;
Smooth = 1;
Hemi = 1;

scrsz = get(0,'ScreenSize');

load dat.mat;
xm = x;

[tfrsq, ConceFT, tfrsqtic] = ConceFT_CWT(time, xm, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi) ;
%[~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
figure;
imageSQ(time, tfrsqtic, abs(ConceFT)); axis square;
xlabel('Time') ; ylabel('Freq') ; 
title('MT=30');
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsCWT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
