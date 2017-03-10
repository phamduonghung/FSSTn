% This code draws the synchrosqueezed wave packet transform in Figure 1 of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

%% set up data
is_real = 0;
load dat.mat;
fff = x.';
N = length(f2);

%% set up parameters
res = 0.2;
NG = N;
is_unif = 1;
typeNUFFT = 1;
is_cos = 1;
epsl = 1e-2;
xo = x;
highglobal = 130;
rad = 1.3;
red = 10;
t_sc = 1/2 + 2/8;
R_low = 0; R_high = highglobal;

%% compute the transform
tic;
T_f = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
toc

%% display
head = sprintf('sswpt energy distribution,red = %d',red);
pic = figure;imagesc([0 1],[0 highglobal],real(T_f((end-1)/2:end,:)));
axis square; axis xy;
T_ftemp = T_f;
xlabel('Time')
% ylabel('Freq')
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsWPT';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);


