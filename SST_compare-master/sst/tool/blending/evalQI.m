function [myQ] = evalQI(tau, f, t, m) ;

M = length(tau) ;
N = M - m - 1 ;
r = N - m - 1 ;

alpha_mjk = evalalpha_mjk(tau, m) ;
myQ = zeros(length(t), 1) ;

for ii = 1: length(t)
	if ~mod(ii, 10000) ; fprintf('*') ; end
	x = t(ii) ;

	Nmjx = zeros(M, 1) ;
	for j = -m : r
		Nmjx(j + m + 1) = eval_Nmj(tau, m, j, x) ;
	end

	NmjxIdx = find(Nmjx ~= 0) ;

		%% get (10.4.34)
	Umkx = zeros(M, 1) ;
	for k = -m: r+m
		if (max(k-m, -m) + m + 1 <= NmjxIdx(end)) & (min(k, r) + m + 1 >= NmjxIdx(1))
			jidx = [ max(k-m, -m) : min(k, r) ] + m + 1 ;
			Umkx(k + m + 1) = alpha_mjk(jidx, k + m + 1)' * Nmjx(jidx) ;
		end
	end
	
	UmkxIdx = find(Umkx ~= 0) ;
		%% get (10.4.33)
	myQ(ii) = f(UmkxIdx) * Umkx(UmkxIdx) ;

end
end


%==================================================
function [alpha_mjk] = evalalpha_mjk(tau, m)
%
%	evaluate alpha_{m,j,k} in (10.4.52)
%
M = length(tau) ;
N = M - m - 1 ;
r = N - m - 1 ;


alpha_mjk = zeros(N, M) ;

for ii = 1: N
	for kk = 0: m
			%% evaluate the denominator	\Pi_{k\neq l=0}^m (\tau_{j+k}-\tau_{j+l})
		den = 1 ; 
		for ll = [0:kk-1 kk+1:m]
			den = den * (tau(ii + kk) - tau(ii + ll)) ;
		end
		
			%% evaluate the numerator
		AllPer = perms([1:m]) ;
		num = 0 ;
		for pp = 1:size(AllPer, 1)
			per = AllPer(pp, :) ;	
			numprod = 1 ;
				%% \Pi_{k\neq l=0}^m (\tau_{j+v_l}-\tau_{j+l})
			vll = 1 ;
			for ll = [0:kk-1 kk+1:m]
				numprod = numprod * (tau(ii + per(vll)) - tau(ii + ll)) ;
				vll = vll + 1 ; 
			end
				%% sum_{{v_0,\ldots,v_m}\{v_k}\subset per{1,\ldots,m}} \Pi_...
			num = num + numprod ;
		end
	
			%% evaluate the number in (10.4.52)
		alpha_mjk(ii, ii+kk) = num/den/factorial(m) ;	
		
	end
end

end



%==================================================
function [Nmjx] = eval_Nmj(tau, m, j, x)
%
% evalute the m-th B spline associated with the knot tau_j at x
%
if length(x) > 1
	error('only support one value at once\n') ;
end

if x < tau(1) | x > tau(end)
	error('Out of the scope\n') ;
end
	
summand = 0 ;
for k = j : j+m+1
	den = 1 ;
	for l = [j:k-1 k+1:j+m+1]
		den = den * (tau(l + m + 1) - tau(k + m + 1)) ;
	end
	summand = summand + ( max(0, x - (tau(k + m + 1))).^m ) / den ;
end

	%% (10.2.32) 
Nmjx = (tau(j + 2*m + 2) - tau(j + m + 1)) * summand ; 

end

