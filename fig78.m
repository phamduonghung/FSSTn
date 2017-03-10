function fig78()
% fig78 : draws Figure 78 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.

clc; clear all; close all;

chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

set(0,'DefaultAxesFontSize',16);

cas=2;
if cas == 1 % Livingston signal
    load fig1-observed-Livingston.txt
    t1 = fig1_observed_Livingston(:,1)';
    s1 = fig1_observed_Livingston(:,2)';
elseif cas ==2 % Hanford signal
    load fig1-observed-Hanford.txt
    t1=fig1_observed_Hanford(:,1)';
    s1=fig1_observed_Hanford(:,2)';
elseif cas ==3 % Hanford theorical waveform
    load fig1-waveform-H.txt
    t1=fig1_waveform_H(:,1)';
    s1=fig1_waveform_H(:,2)';
elseif cas ==4  % Livingston theorical waveform
    load fig1-waveform-L.txt
    t1=fig1_waveform_L(:,1)';
    s1=fig1_waveform_L(:,2)';
else % Hanford  unfilterd waveform
    load fig2-unfiltered-waveform-H.txt
    t1=fig2_unfiltered_waveform_H(:,1)';
    s1=fig2_unfiltered_waveform_H(:,2)';
end

load fig1-waveform-H.txt
load fig1-residual-H.txt
sw=fig1_waveform_H(:,2)';
resi_existing = fig1_residual_H(:,2)';

N1 = length(s1);
Te1=t1(2)-t1(1);

%%zeros pading for signal
nv = log2(N1);
 if mod(nv,1)~=0
     warning('The signal is not a power of two, zero padding to the next power');
     s1=[s1 zeros(1,2^(floor(log2(N1))+1)-N1)];
     t1pad = max(t1)+Te1:Te1:max(t1)+(2^(floor(log2(N1))+1)-N1)*Te1;
     t1=[t1 t1pad];
     sw=[sw zeros(1,2^(floor(log2(N1))+1)-N1)];
 end

N = length(s1);

fs = 0:N/32;

% Parameters
gamma = 10^(-3);
d=10;
sigma=0.05;
index = floor(N1/8+1):floor(7*N1/8);
clwin = 10;
lambda = 0;
nmodes = 1;

% Choice of time and frequency bins
ft =1:N/32;bt=1:N;
 
%% signal
s = s1;
s=hilbert(s);

%% TF presentations of signal

[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigma,ft,bt);

%% Display TFRs
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle1(1) = figure; 
set(FigHandle1(1),'units','normalized','outerposition',[0 0 1 1]);

ha = subplot(2,5,1);
plot(t1(1:N1),real(s(1:N1)),'m','LineWidth',1);
ylabel('Strain (10^{-21})');
xlim([min(t1) 1.01*max(t1(1:N1))]);
ylim([-1 2]);
legend('H1 observed')

FigHandle1(2) = figure; 
set(FigHandle1(2),'units','normalized','outerposition',[0 0 1 1]);

ha = subplot(2,5,2);
imagesc(t1(1:N1),fs,abs(STFT(:,1:N1)));
set(gca,'YDir','normal')
set(ha,'xtick',[],'ytick',[]);xlabel('');ylabel('');

FigHandle1(3) = figure; 
set(FigHandle1(3),'units','normalized','outerposition',[0 0 1 1]);

ha = subplot(2,5,3);
imagesc(t1(1:N1),fs,abs(FSST2(:,1:N1)));%title('(c) FSST2');
set(gca,'YDir','normal')
set(ha,'xtick',[],'ytick',[]);xlabel('');ylabel('');

FigHandle1(4) = figure; %colormap(1-gray);
set(FigHandle1(4),'units','normalized','outerposition',[0 0 1 1]);

ha = subplot(2,5,4);
imagesc(t1(1:N1),fs,abs(FSST4(:,1:N1)));%title('(d) FSST4');
set(gca,'YDir','normal')
set(ha,'xtick',[],'ytick',[]);xlabel('');ylabel('');

%%%%%%%%%%%%%%%%%%%%%% print TFR Hanford
for i =1:4
export_fig(FigHandle1(i), ... % figure handle
    sprintf('%s/TFRs_Hanford_%d', chemin0,i),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
end

%% Ridge extraction FSST2
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle(1) = figure;
set(FigHandle(1), 'Position', [300,1000, 1000, 800]);
h=subplot(2,2,1);
imagesc(t1(1:N1),fs,abs(FSST2(:,1:N1)));
set(gca,'YDir','normal')

set(gca,'ytick',[]);
set(gca,'xtick',[]);

[Cs, Es] = exridge_mult_Noise(FSST2, nmodes, lambda, clwin);
colridge_hanford(h,Cs(:,1:N1),2);
legend('FSST2 ridge extraction','Location','northwest');
xlabel('time'); ylabel('frequency');

imf2 = sigma*recmodes(FSST2,Cs,d);

%% FSST4
FigHandle(2) = figure; 
set(FigHandle(2), 'Position', [300,100, 1000, 800]);
h=subplot(2,2,1);
imagesc(t1(1:N1),fs,abs(FSST4(:,1:N1)));
set(gca,'YDir','normal')
set(gca,'ytick',[]);
set(gca,'xtick',[]);

[Cs, Es] = exridge_mult_Noise(FSST4, nmodes, lambda, clwin);
colridge_hanford(h,Cs(:,1:N1),3);
legend('FSST4 ridge extraction','Location','northwest');
xlabel('time'); ylabel('frequency');

imf4 = sigma*recmodes(FSST4,Cs,d);

%% Upper and Lower Subplots with Titles
 
FigHandle(3) = figure;
set(FigHandle(3), 'Position', [300,100, 1000, 800]);
h=subplot(2,2,1);
plot(t1(1:N1),real(sw(1:N1)),'b',t1(1:N1),real(imf2(1,1:N1)),'r',t1(1:N1),real(imf4(1,1:N1)),'g','LineWidth',1);
xlim([min(t1(1:N1)) 1.01*max(t1(1:N1))]);
ylim([-1 2]);
legend('Numerical relativity','FSST2 reconstructed signal','FSST4 reconstructed signal','Location','northwest');

FigHandle(4) = figure;
set(FigHandle(4), 'Position', [300,100, 1000, 800]);
h=subplot(2,2,1);
plot(t1(1:N1),real(imf2(1,1:N1)-real(sw(1:N1))),'r',t1(1:N1),real(imf4(1,1:N1)-real(sw(1:N1))),'g','LineWidth',1);
xlim([min(t1(1:N1)) 1.01*max(t1(1:N1))]);
ylim([-1 2]);
legend('FSST2 residual','FSST4 residual','Location','northwest');

disp('FSST2 SNR')
disp(snr(real(sw(index)),real(imf2(1,index))-real(sw(index))));
disp(snr(real(resi_existing(index)),real(imf2(1,index)-real(sw(index)))-real(resi_existing(index))));

disp('FSST4 SNR')
disp(snr(real(sw(index)),real(imf4(1,index))-real(sw(index))));
disp(snr(real(resi_existing(index)),real(imf4(1,index)-real(sw(index)))-real(resi_existing(index))));

%%%%%%%%%%%%%%%%%%%%%% print FSST2 Hanford
for i=1:4
export_fig(FigHandle(i), ... % figure handle
    sprintf('%s/reconstruction_Hanford_%d', chemin0,i),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...           % file format
    '-r5000' );           % resolution in dpi
end
end
