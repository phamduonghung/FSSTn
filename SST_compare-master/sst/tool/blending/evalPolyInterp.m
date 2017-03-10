function [interpf] = evalPolyInterp(tau, f, interptau) ;

interpf = zeros(size(interptau)) ;

for ii = 1: length(interptau)
	
	x = interptau(ii) ;
	
	Lnj = zeros(length(tau), 1) ; 
	for j = 0: length(tau)-1
		m = 1 ;
		for k = [0:j-1 j+1:length(tau)-1]
			m = m * (x - tau(k+1)) ./ (tau(j+1) - tau(k+1)) ;
		end
		Lnj(j+1) = m ;
	end
	interpf(ii) = f' * Lnj ;

end

end
