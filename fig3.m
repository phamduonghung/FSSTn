function fig3()

% fig3 : draws Figure 3 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
N = 1024;
fs = 0:N/2;
gamma = 0;
sigma = 0.05;

t = (0:N-1)/N;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
 [a1,a2,if1,if2,s1,s2,s] = signal_test(t);
 
%% TF presentations of signal
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigma,ft,bt);
[~, ~, RM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(s,gamma,sigma);

%% zoom coordinates
xn=0.07;xN=0.18;
yn1=220;yN1=490;
yn=90;yN=120;

%% Display TFRs - Figure 2
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0.05 0.1], [0.05 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle1(1) = figure; %colormap(1-gray);
set(FigHandle1(1),'units','normalized','outerposition',[0 0 1 1]);

subplot(6,6,[1,7,13,19,25,31]);
imagesc(t,fs,abs(STFT)); %title('(a) STFT');
set(gca,'YDir','normal');set(gca, 'xtick', []) ;
rectangle('Position',[xn yn xN-xn yN-yn],'EdgeColor','green','Linewidth',1.5);
rectangle('Position',[xn yn1 xN-xn yN1-yn1],'EdgeColor','red','Linewidth',1.5);

FigHandle1(2) = figure; %colormap(1-gray);
set(FigHandle1(2),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,2);
imagesc(t,fs,abs(STFT));%title('(b) STFT ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(3) = figure; %colormap(1-gray);
set(FigHandle1(3),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,8);
imagesc(t,fs,abs(RM));%title('(c) RM ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(4) = figure; %colormap(1-gray);
set(FigHandle1(4),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,14);
imagesc(t,fs,abs(FSST));%title('(d) FSST ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(5) = figure; %colormap(1-gray);
set(FigHandle1(5),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,20);
imagesc(t,fs,abs(FSST2));%title('(e) FSST2 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(6) = figure; %colormap(1-gray);
set(FigHandle1(6),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,26);
imagesc(t,fs,abs(FSST3));%title('(f) FSST3 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(7) = figure; %colormap(1-gray);
set(FigHandle1(7),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,32);
imagesc(t,fs,abs(FSST4));%title('(g) FSST4 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn yN]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'green';
ha.XColor = 'green';

FigHandle1(8) = figure; %colormap(1-gray);
set(FigHandle1(8),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,3);
imagesc(t,fs,abs(STFT)); %title('(h) STFT ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

FigHandle1(9) = figure; %colormap(1-gray);
set(FigHandle1(9),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,9);
imagesc(t,fs,abs(RM)); %title('(k) RM ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

FigHandle1(10) = figure; %colormap(1-gray);
set(FigHandle1(10),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,15);
imagesc(t,fs,abs(FSST)); %title('(l) FSST ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

FigHandle1(11) = figure; %colormap(1-gray);
set(FigHandle1(11),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,21);
imagesc(t,fs,abs(FSST2)); %title('(m) FSST2 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

FigHandle1(12) = figure; %colormap(1-gray);
set(FigHandle1(12),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,27);
imagesc(t,fs,abs(FSST3)); %title('(n) FSST3 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

FigHandle1(13) = figure; %colormap(1-gray);
set(FigHandle1(13),'units','normalized','outerposition',[0 0 1 1]);

ha=subplot(6,6,33);
imagesc(t,fs,abs(FSST4)); %title('(o) FSST4 ');
set(gca,'YDir','normal')
set(ha,'xlim',[xn xN],'ylim',[yn1 yN1]);
set(ha,'xtick',[],'ytick',[],'Linewidth',1.5);xlabel('');ylabel('');
ha.YColor = 'red';
ha.XColor = 'red';

%%%%%%%%%%%%%%%%%%%%%% print TFR Figure 2

for i = 1:13
    export_fig(FigHandle1(i), ... % figure handle
        sprintf('%s/TFRs_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end
end
