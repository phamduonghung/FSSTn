	%% test the non-uniform sampling blending spline
clear; close all;

n0 = 2048 ;
	%% spline order. m = 3 means cubic spline
m = 3 ;
sr = 64 ;

	% truth
%t = linspace(1/n0, n0/sr, n0) ;
t = linspace(0, (n0+1)/sr, n0+2) ;
x = randn(n0+2, 1) ;
x = smooth(x, 512, 'loess') ;

	%% sampling
idx = [1 sr/2:sr:n0 n0] ;
P = 48 ;
perp = floor(P * ( rand(1, length(idx)-2) - .5 )) ;  
idx(2:end-1) = idx(2:end-1) + perp ;
tau0 = t(idx) ;
f0 = x(idx) ;
r = length(f0) ;

gap = median(tau0(2:end) - tau0(1:end-1)) ;
	%% tau0(1) should be normalized to 1/sr 
tau = [ [-m:-1]*gap tau0 tau0(end) + [1:m]*gap ] ;
[L] = evalPolyInterp(tau0(1:m+1), f0(1:m+1), [-m:-1]*gap) ;
[R] = evalPolyInterp(tau0(r-m:r), f0(r-m:r), tau0(end) + [1:m]*gap) ; 
%f = [ [-m:-1]*0 0 f0' 0 [1:m]*0 ] ; 
f = [ L f0' R ] ; 

	%% my Quasi-Interpolation  
fprintf('Get QI...') ; t1 = tic ; 
myQ = evalQI(tau, f, t, m) ;
fprintf([num2str(toc(t1)),'\n']) ;

	%% my R operator
fprintf('Get RI...') ; t2 = tic ; 
Qf0 = myQ(idx) ;
[L] = evalPolyInterp(tau0(1:m+1), Qf0(1:m+1), [-m:-1]*gap) ;
[R] = evalPolyInterp(tau0(r-m:r), Qf0(r-m:r), tau0(end) + [1:m]*gap) ;
%Qf = [ [-m:-1]*0 0 Qf0' 0 [1:m]*0 ] ;
Qf = [ L Qf0' R ] ;
myR = evalRI(tau, f - Qf, t, m) ;  
fprintf([num2str(toc(t2)),'\n']) ;

	%% my blending operator
myB = myQ + myR  ;

	%% cubic spline provided by matlab
matlabF = interp1(tau0, f0, t, 'spline', 'extrap') ;


	%% plot the results
plot(t(:), x(:), 'color', [.8 .8 .8], 'linewidth', 1.2) ; 
hold on ;
plot(tau0(:), f0(:), 'ko', 'linewidth', 2) ;
plot(t(:), matlabF(:), '-.', 'color', [.0 .0 .0]) ; 
plot(t(:), myB(:), 'r--', 'linewidth', 1.2) ;
set(gca, 'fontsize', 18) ; 
legend('truth', 'sampling', 'matlab cubic', 'blending');
axis tight;
