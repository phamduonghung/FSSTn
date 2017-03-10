function [gammal] = evalBSplineCoeff(f, m, k)
	%% now only even m and k=1 are supported
	%% cubic spline is m=4

if mod(m, 2)
	error('Currently only even m is supported\n') ;
end

if k ~= 1
	error('Currently only k=1 is supported\n') ;
end

vn = [1/48 -1/12 -1/8 7/12 29/24 7/12 -1/8 -1/12 1/48] ;  % (4.6.28)
gammal = conv(upsampling(f), vn, 'same') ;	% (4.6.30)	

end

function [ajj] = upsampling(aj)
    ajj = zeros(length(aj)*2, 1);
    ajj(2:2:end) = aj ;
end
