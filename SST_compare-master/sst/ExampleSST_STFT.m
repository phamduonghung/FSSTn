clear ; close all ;

addpath('./tool') ;

scrsz = get(groot,'ScreenSize');

	%% generate the simulated data
	%% she sampling time (100Hz sampling rate)
	%% high sampling rate to avoid sampling issue
    %% You can generate other different kinds using genSignal.m

Hz = 100 ;
L = 12 ;
N = Hz * L ;


opts.WinLen = 301 ;
        %% Setup parameters
        %% alpha is the grid size in the freq axis
alpha = 0.01/Hz ;


initstate(1) ;


	%% the amplitude modulation of the simulated signal
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;
am1(1:300) = 0 ;
am2(end-400:end) = 0 ;


    %% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = 10 + 6 * if1 ./ max(abs(if1)) ;
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
if2 = pi + 3 * if2 ./ max(abs(if2)) ;
phi1 = cumsum(if1) / Hz ; 
phi2 = cumsum(if2) / Hz ; 
if1(1:300) = nan ;
if2(end-400:end) = nan ;

	%% the simulated signal.
clean = am1 .* cos(2*pi*phi1) + am2 .* cos(2*pi*phi2) ; 
	%% add noise (Gaussian white noise)

noise = 0 * randn(N, 1) ; 
%dis = random('t',4,N, 1) ;
%e = armaxfilter_simulate(dis, .5, 1, .5, 1, -.5) ;
%noise = sigma * e ./ std(e) ;
snrdb = 20 * log10(std(clean)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;


	%% simulated observed time series
xm = clean + noise ;


time = [1/Hz:1/Hz:L]' ;


%%%% call sqSTFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h,Dh,t] = hermf(opts.WinLen, 1, 6) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(xm, 0/Hz, 20/Hz, 0.5/length(xm)/2, 1, h', Dh') ; 

subplot(1,2,1) ;
imageSQ(time, tfrsqtic*Hz, log(1+abs(tfrsq))) ; colormap(1-gray) ;
set(gca,'fontsize', 20) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
title('SST on STFT\newline(displayed in the log scale)') ;

subplot(1,2,2) ;
imageSQ(time, tfrsqtic*Hz, abs(tfrsq)) ; colormap(1-gray) ;
set(gca,'fontsize', 20) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
title('SST on STFT\newline(displayed in the linear scale)') ;

