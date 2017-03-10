function [t, s, s1, s2, am1, am2, if1, if2, cf1, cf2, trend, noise] = genSignal(type, sigma) ;


if type == 1

    N = 512; % length of signal in sample
    Fs = 50;
	t = [1/Fs:1/Fs:N/Fs]' ;
	
	am1 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess') ;
	am2 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess') ;
	am1 = 1 + (am1 + 2*max(abs(am1))) ./ (3*max(abs(am1))) ;
	am2 = 1 + (am2 + 2*max(abs(am2))) ./ (3*max(abs(am2))) ;
	
		% the instantaneous frequency of the simulated signal
	if1c = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess') ;
	if2c = smooth(cumsum(randn(N,1)) ./ Fs, 100, 'loess') ;
	if1 = 6 + 6*(if1c + 1.2*max(abs(if1c))) ./ ( 2.2*max(abs(if1c)) ) + t ;
	if2 = pi + 6*(if2c + 1.2*max(abs(if2c))) ./ ( 2.2*max(abs(if2c)) ) ;
	cf1 = [(if1(2:end)-if1(1:end-1)) .* Fs] ; cf1 = [cf1; cf1(end)] ;
	cf2 = [(if2(2:end)-if2(1:end-1)) .* Fs] ; cf2 = [cf2; cf2(end)] ;
	phi1 = cumsum(if1) / Fs ;
	phi2 = cumsum(if2) / Fs ;
	
	mask1 = ones(N, 1) ; mask1(1:102) = 0 ;
	mask2 = ones(N, 1) ; mask2(end-102:end) = 0 ;
	
		% the simulated signal.
	s1 = am1 .* cos(2*pi* phi1) .* mask1 ;
	s2 = am2 .* cos(2*pi* phi2) .* mask2 ;
	clean = s1 + s2 ;
	if1(1:102) = nan ; if2(end-102:end) = nan ;
	am1(1:102) = nan ; am2(end-102:end) = nan ;

	trend = 0 * (4 + smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess')) ;
	
		% add noise (Gaussian white noise) 2 = -3dB, 1= 3dB
	noise = sigma * randn(N, 1) ;
	snrdb = 20 * log10(std(clean) ./ std(noise)) ;
	
		% simulated observed time series
	s = clean + trend + noise ;

	fprintf(['Singal type ',num2str(type),'; snrdb = ',num2str(snrdb),'\n']) ;

%===============================================
elseif type == 2
	
    N = 512; % length of signal in sample
    Fs = 50;
    t = [1/Fs:1/Fs:N/Fs]' ;

    am1 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess') ;
    am2 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess') ;
    am1 = 1 + (am1 + 2*max(abs(am1))) ./ (3*max(abs(am1))) ;
    am2 = 1 + (am2 + 2*max(abs(am2))) ./ (3*max(abs(am2))) ;

		%% 2 linear chirps
	mask1 = ones(N, 1) ; mask1(1:100) = 0 ;
	mask2 = ones(N, 1) ; 
	s1 = am1 .* cos(2*pi*(3*t + t.^2)) .* mask1 ;
	s2 = am2 .* cos(2*pi*(20*t - t.^(1.8))) ;
	if1 = (3*ones(N,1) + 2*t) .* mask1 ;
	if2 = 20*ones(N,1) - 1.8*t.^(0.8) ;
	cf1 = 2*ones(N,1) .* mask1 ;
	cf2 = 0.8*1.8*t.^(-0.2) ;
	s = s1 + s2 ;
	noise = sigma * randn(size(s)) ; trend = zeros(size(s)) ;
	am1(1:100) = nan ; if1(1:100) = nan ;
	snrdb = 20 * log10(std(clean) ./ std(noise)) ;

	fprintf(['Singal type ',num2str(type),'; snrdb = ',num2str(snrdb),'\n']) ;


%===============================================
elseif type == 3

		% Semi Real Signal
	load('ExampleAf.mat');
	time = t(401:100:end-100) ;
	Fs = 10 ;
	
	am1 = smooth(cumsum(randn(length(time),1)) ./ 32, 300, 'loess') ;
	am1 = 1 + (am1 + max(abs(am1))) ./ (3*max(abs(am1))) ;
	phi1 = cumsum(rri(401:100:end-100))/Fs ;
	if1 = resample(rri(401:100:end-100),1,1);
	cf1 = (if1(2:end) - if1(1:end-1)) * Fs ;
	cf1 = [cf1; cf1(end)] ;
	
	am2 = smooth(cumsum(randn(length(am1),1)) ./ Fs, 300, 'loess') ;
	am2 = 1 + (am2 + max(abs(am2))) ./ (3*max(abs(am2))) ;
	if2c = smooth(cumsum(randn(length(am1),1)) ./ Fs, 20, 'loess') ;
	if2 = 2.5 + 2*(if2c + 0.5*max(abs(if2c))) ./ ( 1*max(abs(if2c)) ) ;
	phi2 = cumsum(if2) / Fs ;
    cf2 = (if2(2:end) - if2(1:end-1)) * Fs ;
    cf2 = [cf2; cf2(end)] ;
	
	mask2 = ones(length(am1), 1) ; mask2(1:202) = 0 ;
	if2 = if2 - sin(time) ;
	if2(1:202) = nan ;
	
	s1 = am1 .* cos(2*pi*phi1) ;
	s2 = am2 .* cos(2*pi*(phi2+cos(time))) .* mask2 ;
	clean = s1 + s2 ;
	
	        % add noise (Gaussian white noise)
	noise = sigma * randn(length(clean), 1) ;
	snrdb = 20 * log10(std(clean) ./ std(noise)) ;
	trend = zeros(size(noise)) ;
	s = clean + noise + trend ;
	t = time - time(1) ;
    N = length(s) ;

	fprintf(['Singal type ',num2str(type),'; snrdb = ',num2str(snrdb),'\n']) ;

%===============================================
elseif type == 4

    	%% try the bad signal
	load badresp
	t = t(1:2:end)' ; s = x(1:2:end) ;
	Fs = 1./(t(2)-t(1)) ;
	N = length(s) ;

	s1 = zeros(size(t)) ; s2 = zeros(size(t)) ;
	am1 = zeros(size(t)) ; am2 = zeros(size(t)) ;
	if1 = zeros(size(t)) ; if2 = zeros(size(t)) ;
	cf1 = zeros(size(t)) ; cf2 = zeros(size(t)) ;
	noise = zeros(size(t)) ; trend = zeros(size(t)) ;

	snrdb = 20 * log10(std(clean) ./ std(noise)) ;
	fprintf(['Singal type ',num2str(type),'; snrdb = ',num2str(snrdb),'\n']) ;

end

