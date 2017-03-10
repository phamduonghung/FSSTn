function [myR] = evalRI(tau, f, t, m) ;

if m ~= 3
	error('Now only support m = 3\n') ;
end 

M = length(tau) ;
N = M - m - 1 ;
r = N - m - 1 ;

myR = zeros(length(t), 1) ;

M2 = M*2 - 1;
tau2 = zeros(M2, 1) ;
for ii = 1:length(tau) - 1
	tau2((ii-1)*2 + 1) = tau(ii) ;
	tau2(ii*2) = (tau(ii) + tau(ii+1)) / 2 ;   
end
tau2(end) = tau(end) ;

	%% only hold for the cubic blending
Delta = 2 ;

for ii = 1: length(t)
	if ~mod(ii, 10000) ; fprintf('*') ; end

	x = t(ii) ;

	Nmjx = zeros(M, 1) ;
	for j = 1 : r
		jidx = 2*m + 2*j + 1 - Delta : 2*m + 2*j + 1 + Delta ;
		tauS = tau2(2*m + 2*j + 1 - Delta) ;
		posden = tau2(2*m + 2*j + 1) ;
		Nmjx(j + m + 1) = eval_Nmj2(tau2(jidx), x) ./ eval_Nmj2(tau2(jidx), posden) ;
	end

%if ~mod(ii, 256); keyboard; end
	NmjxIdx = find(Nmjx ~= 0) ;
		%% get (10.4.33)
	myR(ii) = f(NmjxIdx) * Nmjx(NmjxIdx) ;

end
end



%==================================================
function [Nmjx] = eval_Nmj2(tau, x)
%
% evalute the m-th B spline associated with the knot tau_j at x
%
if length(x) > 1
	error('only support one value at once') ;
end

m = 3 ;

	
summand = 0 ;
for k = 1 : length(tau)
	den = 1 ;
	for l = [1:k-1 k+1:length(tau)]
		den = den * (tau(l) - tau(k)) ;
	end
	summand = summand + ( max(0, x - (tau(k))).^m ) / den ;
end

	%% (10.2.32) 
Nmjx = (tau(end) - tau(1)) * summand ; 

end

