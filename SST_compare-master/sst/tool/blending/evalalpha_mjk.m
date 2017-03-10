function [alpha_mjk] = evalalpha_mjk(tau, m)
	%% M = | [-m:1:r+1+m] | ;
M = length(tau) ;
N = M - m - 1 ;
r = N - m - 1 ;


alpha_mjk = zeros(M, M) ;

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


