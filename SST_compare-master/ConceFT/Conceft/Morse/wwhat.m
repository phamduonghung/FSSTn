function vector=wwhat(f,beta,gam,k)
% This function calculates morse complex wavelets (in frequency domain)
% It is very similar to morsefreqs
c=-1+((2.*beta+1)./gam);
con=sqrt(gam.*(2.^((2.*beta+1)./gam)).*gamma(k+1));
con2=sqrt(pi.*gamma(k+((2.*beta+1)./gam)));
con=sqrt(pi./2).*(con./con2);
arr=2.*((2.*pi.*abs(f)).^gam);
m=con.*(2.*pi.*abs(f)).^(beta).*exp(-(2.*pi.*abs(f)).^gam).*laggen(arr,k,c);
m=sqrt(4.*pi).*m;
for j=1:length(f)
if f(j)<0
m(j)=0;
end;
end;
vector=m;