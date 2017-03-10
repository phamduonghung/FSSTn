function [wmk] = evalwmk(m)
	%% evaluate w_{m,k} = N_m(k) in (4.3.10) using (4.1.12)

k = zeros(m+1);
m1 = zeros(1, m+1);
for ii = 1:m+1
	m1(ii) = ((-1)^(ii-1)) * choose(m, ii-1) ;
end

for ii = 1:m+1
	k(ii, ii:end) = [0:m+1-ii] .^ (m-1) ;
	k(ii, :) = m1 .* k(ii, :) ;
end

wmk = sum(k, 2) ./ factorial(m-1) ;
wmk = wmk( 2: end-1 );

end

