function [c] = evalSplineCoeff(f, m)
	%% now only even m is supported
	%% cubic spline is m=4

if mod(m, 2)
	error('Currently only even m is supported\n');
end

J = length(f) ;
wmk = evalwmk(m) ;
N = zeros(J) ;

for ii = 1: m/2
	N = N + wmk(m/2 + ii - 1) * diag(ones(J - ii + 1, 1), ii - 1) ;
	if ii > 1
		N = N + wmk(m/2 + ii - 1) * diag(ones(J - ii + 1, 1), 1 - ii) ;
	end
end

c = inv(N) * f(:) ;

end
