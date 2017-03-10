function fig1()
% fig1 : draws Figure 1 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.
%
clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
N = 1024;

t = (0:N-1)/N;

%% Test signal 
 [a1,a2,if1,if2,s1,s2,s] = signal_test(t);
 
%% display waveform signal - Figure 1 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle(1) = figure; 
set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s1), 'm'); 
hold on; plot(t, a1, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-10 10]);
legend('Re(f_1)','A_{1}','Location','south','Orientation','horizontal') ; 

FigHandle(2) = figure; 
set(FigHandle(2),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s2), 'g'); 
hold on; plot(t, a2, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-10 10]);
legend('Re(f_2)','A_{2}','Location','south','Orientation','horizontal') ; 

FigHandle(3) = figure; 
set(FigHandle(3),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s), 'b') ;
legend('Re(f)','Location','south') ; 
set(ha,'ylim',[-20 20]);

%%%%%%%%%%%%%%%%%%%%%% print Figure 1

for i=1:3
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/MCS_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end
end
