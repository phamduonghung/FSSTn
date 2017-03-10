function [F] = evalSplineInterp(c, m, J)
	%% increase the resoluation from 1 to 2^J 

	pmk = evalpmk(m) ;
	wmk = evalwmk(m) ;

	aj = c ; 
	for jj = 1:J
		ajj = upsampling(aj) ;			% (4.3.11)
		aj = conv(ajj, pmk, 'same') ;	% (4.3.12)
	end

	F = conv(aj, wmk, 'same') ;

end

function [ajj] = upsampling(aj)
	ajj = zeros(length(aj)*2, 1);
	ajj(2:2:end) = aj ;
end 
