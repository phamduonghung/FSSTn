function [pmk] = evalpmk(m)

pmk = zeros(m+1, 1);
for ii = 1:m+1
	pmk(ii) = ( 2^(-m+1) ) * choose(m, ii-1) ;
end

end


