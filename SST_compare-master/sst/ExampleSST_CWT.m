clear ; %close all ;

addpath('./tool') ;


%% fix the seed
initstate(23400) ;
opts = struct();
opts.motherwavelet = 'Cinfc' ;
opts.CENTER = 1 ;
opts.FWHM = 0.3 ;

lowfreq = 0 ;
highfreq = 10 ;

%% high sampling rate to avoid sampling issue
Hz = 50 ;
L = 20 ;
time = [1/Hz:1/Hz:L]' ;
N = length(time) ;
alpha = 0.05 ; 

scrsz = get(0,'ScreenSize');



%=====================================================
	%% generate simulated signals. You can generate more using genSignal.m

%% the amplitude modulation of the simulated signal
%% if you don't have statistics toolbox, try to use convolution to do smoothing
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;
am1(1:200) = 0 ;
am2(end-230:end) = 0 ;

%% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = 5 + if1 ./ max(abs(if1)) ;
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
if2 = pi + if2 ./ max(abs(if2)) ;
phi1 = cumsum(if1) / Hz ;
phi2 = cumsum(if2) / Hz ;

%% the simulated signal.
clean = am2 .* cos(2*pi*(phi2)) + am1 .* cos(2*pi*phi1) ; %+ 2*cos(2*pi*(time+0.4*time.^2));
if2 = if2 ;% - 0.2*time ;
if1(1:200) = nan ;
if2(end-230:end) = nan ;
if3 = 1+0.8*time ;

%% add noise (Gaussian white noise)
noise = 2 * randn(N, 1) ;
noise(1:end/2) = 1*randn(N/2,1) ;
xm = clean ; % + noise ;


%==================

[~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);




%%
figure

subplot(1,2,1);
imageSQ(time, tfrsqtic, log(1+abs(tfrsq))); 
xlabel('Time (sec)') ; ylabel('Freq') ; title('SST (display in log scale)') ;
subplot(1,2,2);
imageSQ(time, tfrsqtic, abs(tfrsq)); 
xlabel('Time (sec)') ; ylabel('Freq') ; title('SST (display in linear scale)') ;
