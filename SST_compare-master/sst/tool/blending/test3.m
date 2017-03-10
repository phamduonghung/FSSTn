	%% test the non-uniform sampling blending spline
clear; close all;

n0 = 2048 ;
	%% spline order. m = 3 means cubic spline
m = 3 ;
sr = 64 ;

	% truth
t = linspace(0, (n0+1)/sr, n0+2) ;
x = randn(n0+2, 1) ;
x = smooth(x, 512, 'loess') ;

	%% sampling
idx = [1 sr/2:sr:n0 n0] ;
P = ceil(sr*2/3) ;
perp = floor(P * ( rand(1, length(idx)-2) - .5 )) ;  
idx(2:end-1) = idx(2:end-1) + perp ;
	%% tau0(1) should be normalized to 0 
tau0 = t(idx) ;
f0 = x(idx) ;
    %% cubic spline provided by matlab
matlabF = interp1(tau0, f0, t, 'spline', 'extrap') ;


%===================================================
	%% Simulation 1: get all interpolations at once (not real time)
[myB] = evalBlending(tau0, f0, t, m) ;

	%% plot the results
plot(t(:), x(:), 'color', [.8 .8 .8], 'linewidth', 1) ; 
hold on ;
plot(tau0(:), f0(:), 'ko', 'linewidth', 1) ;
plot(t(:), matlabF(:), '-.', 'color', [.0 .0 .0]) ; 
plot(t(:), myB(:), 'r--', 'linewidth', 1) ;
set(gca, 'fontsize', 18) ; 
legend('truth', 'sampling', 'matlab cubic', 'blending');
axis tight;



%===============================================================
	%% Simulation 2: real time interpolation
	
	%% we only interpolate the data with related to the most recent InterpLen observations
InterpLen = 2*m ;

for rtii = InterpLen+1 : length(tau0)
	tau0rt = tau0(rtii-InterpLen: rtii) - tau0(rtii-InterpLen) ;
	f0rt = f0(rtii-InterpLen: rtii) ;
	trt = t(idx(rtii-InterpLen): idx(rtii)) - t(idx(rtii-InterpLen)) ; 
	[myBrt] = evalBlending(tau0rt, f0rt, trt, m) ;
	matlabFrt = interp1(tau0rt, f0rt, trt, 'spline', 'extrap') ;

	plot(t(idx(rtii-1): idx(rtii)), myBrt(end-(idx(rtii)-idx(rtii-1)):end), 'r', 'linewidth', 1.2) ;
	axis tight; pause
end
