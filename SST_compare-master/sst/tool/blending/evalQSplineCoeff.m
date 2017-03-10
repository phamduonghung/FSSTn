function [c] = evalQSplineCoeff(f, m, k)
	%% now only even m and k=1 are supported
	%% cubic spline is m=4

if mod(m, 2)
	error('Currently only even m is supported\n') ;
end

if k ~= 1
	error('Currently only k=1 is supported\n') ;
end

lambdak = [-1/6 8/6 -1/6]' ;
c = conv(f, lambdak, 'same') ;

end
