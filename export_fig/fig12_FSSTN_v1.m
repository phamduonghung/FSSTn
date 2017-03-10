function fig12_FSSTN()

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
%chemin0 = '~/Dropbox/Notetaking/Journal_IEEE/Paper_Pham_Meignen/Figures';
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
N = 1024;
fs = 0:N/2;
gamma = 0;
sigma = 0.05;
%index = N/8+1:7*N/8;
%d = 1;
%clwin = 10;
%lambda = 0;
%nmodes = 2;

t = (0:N-1)/N;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
 [a1,a2,if1,if2,s1,s2,s] = signal_test(t);
 
%% display waveform signal - Figure 1 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle111 = figure; 
set(FigHandle111,'units','normalized','outerposition',[0 0 1 1]);

ha = subplot(4,3,1) ;
plot(t, real(s1), 'm'); 
hold on; plot(t, a1, 'k', 'linewidth', 2) ;
set(gca, 'xtick', []) ;
set(ha,'ylim',[-10 10]);
legend('Re(f_1)','A_{1}','Location','south','Orientation','horizontal') ; 
title('(a) Mode 1','fontsize', 14); 
set(gca, 'xtick', []) ;

ha = subplot(4,3,4) ;
plot(t, real(s2), 'g'); 
hold on; plot(t, a2, 'k', 'linewidth', 2) ;
set(gca, 'xtick', []) ; 
set(ha,'ylim',[-10 10]);
legend('Re(f_2)','A_{2}','Location','south','Orientation','horizontal') ; 
title('(b) Mode 2','fontsize', 14); 
set(gca, 'xtick', []) ;


ha = subplot(4,3,7) ;
plot(t, real(s), 'b') ;
legend('Re(f)','Location','south') ; 
title('(c) MCS f','fontsize', 14)
set(gca, 'xtick', []) ;
set(ha,'ylim',[-20 20]);

% ha = subplot(4,3, 10) ;
% plot(t, if1, 'm', 'linewidth', 2) ;  hold on ;
% plot(t, if2, 'g', 'linewidth', 2) ;
% text(0.4, if1(200)+30, '$$\phi''_1$$', 'Interpreter','latex', 'fontsize', 16) ;
% text(0.45, if2(500)+100, '$$\phi''_2$$', 'Interpreter','latex', 'fontsize', 16);
% set(ha,'ylim',[0 N/2]);
% title('(d) Ideal IFs','fontsize', 14)

%%%%%%%%%%%%%%%%%%%%%% print Figure 1
%set(FigHandle1, 'Color', 'white'); % white bckgr
%saveas(FigHandle111, sprintf('%s/MCS_1.pdf', chemin0), 'pdf');

export_fig(FigHandle111, ... % figure handle
    sprintf('%s/MCS', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi


%% TF presentations of signal
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(s,gamma,sigma,ft,bt);
[~, ~, RM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(s,gamma,sigma);

%% zoom coordinates
xn=0.07;xN=0.18;
yn1=220;yN1=490;
yn=90;yN=120;

%% Display TFRs - Figure 2
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0.05 0.1], [0.05 0.01]);
if ~make_it_tight,  clear subplot;  end
FigHandle1 = figure; %colormap(1-gray);
%set(FigHandle1,'units','normalized','outerposition',[0 0 1 1]);
set(FigHandle1, 'Position', [300,300, 850, 850]);

subplot(6,3,[1,4,7,10,13,16]);
imagesc(t,fs,abs(STFT)); title('(a) STFT');
set(gca,'YDir','normal')
rectangle('Position',[xn yn xN-xn yN-yn],'EdgeColor','blue','Linewidth',1.5);
rectangle('Position',[xn yn1 xN-xn yN1-yn1],'EdgeColor','red','Linewidth',1.5);

ha=subplot(6,3,2);
imagesc(t,fs,abs(STFT));title('(b) STFT ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';

ha=subplot(6,3,5);
imagesc(t,fs,abs(RM));title('(c) RM ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';

ha=subplot(6,3,8);
imagesc(t,fs,abs(FSST));title('(d) FSST ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';


ha=subplot(6,3,11);
imagesc(t,fs,abs(FSST2));title('(e) FSST2 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';

ha=subplot(6,3,14);
imagesc(t,fs,abs(FSST3));title('(f) FSST3 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';


ha=subplot(6,3,17);
imagesc(t,fs,abs(FSST4));title('(g) FSST4 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'blue';
ha.XColor = 'blue';


ha=subplot(6,3,3);
imagesc(t,fs,abs(STFT)); title('(h) STFT ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

ha=subplot(6,3,6);
imagesc(t,fs,abs(RM)); title('(k) RM ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

ha=subplot(6,3,9);
imagesc(t,fs,abs(FSST)); title('(l) FSST ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';


ha=subplot(6,3,12);
imagesc(t,fs,abs(FSST2)); title('(m) FSST2 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

ha=subplot(6,3,15);
imagesc(t,fs,abs(FSST3)); title('(n) FSST3 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

ha=subplot(6,3,18);
imagesc(t,fs,abs(FSST4)); title('(o) FSST4 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';


%%%%%%%%%%%%%%%%%%%%%% print TFR Figure 2
%set(FigHandle1, 'Color', 'white'); % white bckgr
export_fig(FigHandle1, ... % figure handle
    sprintf('%s/TFRs', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
end
