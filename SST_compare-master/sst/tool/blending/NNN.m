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

