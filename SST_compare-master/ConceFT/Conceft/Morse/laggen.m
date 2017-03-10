function fn=laggen(x,k,c)
%
% compute generalized Laguerre poly L_k^c(x)
%
%
l=length(x);
sn=-1;
s(1:l)=zeros;
for m=0:k
	sn=-sn;
	ga=gammaln(k+c+1)-gammaln(k-m+1)-gammaln(c+m+1)-gammaln(m+1);
	ga=exp(ga);
	s=s+sn.*ga.*x.^m;
end;
fn=s;
	
	
	
