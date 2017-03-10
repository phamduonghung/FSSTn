	%% test the uniform sampling blending spline
clear; close all;

n0 = 2048 ;
	%% display resolution
J = 5 ;
	%% spline order. m = 4 means cubic spline
m = 4 ;
sr = 128 ;
ds = 128 ;

	% truth
t0 = linspace(1/n0, n0/sr, n0) ;
x = randn(n0, 1) ;
x = smooth(x, 512, 'loess') ;
%x = (linspace(-1, 1, n0)).^2 ;
%x = cos(2*pi*2*t0) ;

	%% sampling
idx = ds:ds:n0 ;
P = 0 ;
perp = floor(P * ( rand(1, length(idx)-1) - .5 )) ;  
idx(1:end-1) = idx(1:end-1) + perp ;
f = x(idx) ;
N = length(f) ;

	%% my cubic spline 
[c] = evalSplineCoeff(f, m) ;
[myF] = evalSplineInterp(c, m, J) ;
tF = linspace(0, N, N* 2^J) ; 

	%% blending operator
[c] = evalBSplineCoeff(f, m, 1) ;
tB = linspace(0, N, N* 2^(J+1)) ; 
[myB] = evalSplineInterp(c, m, J) ;

	%% cubic spline provided by matlab
matlabF = interp1(t0(idx), f, t0, 'spline', 'extrap') ;


	%% plot the results
plot(t0(:), x(:), 'color', [.8 .8 .8], 'linewidth', 1.2) ; 
hold on ;
plot(t0(:), matlabF(:), '-.', 'color', [.0 .0 .0]) ; 
plot(tF(:), myF(:), 'linewidth', 1.2) ;
plot(tB(:), myB(:), 'r--', 'linewidth', 1.2) ;
plot(t0(idx), f, 'ko', 'linewidth', 2) ;
set(gca, 'fontsize', 18) ; 
legend('truth', 'matlab cubic', 'my cubic', 'blending', 'sampling');
axis tight;
