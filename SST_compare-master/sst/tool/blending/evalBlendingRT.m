function [RTmyB] = evalBlendingRT(tau0, f0, t2, m)

if m ~= 3
	error(['Only m = 3 is supported now']) ;
end

%===================================================
RTmyB = zeros(size(t2)) ; 
t2idx = find( (t2 <= tau0(2*m+1)) & (t2 >= tau0(1)) ) ; 
if length(t2idx)
	[tmp] = evalBlending(tau0(1:2*m+1), f0(1:2*m+1), t2(t2idx), m) ;
	RTmyB(t2idx) = tmp ;
end

for jj = 2*m+2 : length(tau0)
	t2idx = find( (t2 <= tau0(jj)) & (t2 >= tau0(jj-2*m)) ) ; 
	if length(t2idx)
		[tmp] = evalBlending(tau0(jj-2*m:jj), f0(jj-2*m:jj), t2(t2idx), m) ;
		t2idx2 = find( (t2 <= tau0(jj)) & (t2 >= tau0(jj-m)) ) ;
		RTmyB(t2idx) = tmp(end-length(t2idx)+1:end) ;
	end
end

