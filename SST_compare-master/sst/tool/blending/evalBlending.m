function [myB] = evalBlending(tau0, f0, t, m)
%
% Given data (tau0, f0), this code evaluate the interpolation on the grid t
%

if t(1) ~= 0 
	tau0 = tau0 - t(1) ;
	t = t - t(1) ; 
end

r = length(f0) ;

t2 = union(t, tau0) ;
idx = zeros(r, 1) ;
for ii = 1: r
	idx(ii) = find(t2 == tau0(ii)) ;
end

gap = median(tau0(2:end) - tau0(1:end-1)) ;
tau = [ [-m:-1]*gap tau0 tau0(end) + [1:m]*gap ] ;
[L] = evalPolyInterp(tau0(1:m+1), f0(1:m+1), [-m:-1]*gap) ;
[R] = evalPolyInterp(tau0(r-m:r), f0(r-m:r), tau0(end) + [1:m]*gap) ;
f = [ L f0' R ] ;

    %% my Quasi-Interpolation  
%fprintf('Get QI...') ; t1 = tic ;
myQ = evalQI(tau, f, t2, m) ;
%fprintf([num2str(toc(t1)),'\n']) ;

    %% my R operator
%fprintf('Get RI...') ; ttic = tic ;
Qf0 = myQ(idx) ;
[L] = evalPolyInterp(tau0(1:m+1), Qf0(1:m+1), [-m:-1]*gap) ;
[R] = evalPolyInterp(tau0(r-m:r), Qf0(r-m:r), tau0(end) + [1:m]*gap) ;
Qf = [ L Qf0' R ] ;
myR = evalRI(tau, f - Qf, t, m) ;
%fprintf([num2str(toc(ttic)),'\n']) ;

ttt = zeros(size(t));
for ii = 1:length(t) ; ttt(ii) = find(t2 == t(ii)) ; end
myB = myQ(ttt) + myR ;
    %% the interpolation by the blending operator
%myB = myQ + myR  ;


end
