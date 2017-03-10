function vector=analyt6(beta,gam,Area,nk)

n=2000;

A=(1:100)./15;
j=45;


r=(2*beta+1)/gam;
C=1+Area*gam^2*exp(2*gammaln(r)-gammaln(r+(1/gam))-gammaln(r-(1/gam)))/(2*pi*(r*gam-1));
for k=0:nk-1
lambda(k+1)=betainc((C-1)/(C+1),k+1,r-1)^2;
psihat(1:n/2)=sqrt(A(j)).*wwhat(A(j).*[(0:n/2-1)/(n)],beta,gam,k);
psihat(n/2+1:n)=sqrt(A(j)).*wwhat(A(j).*[(-n/2:-1)/(n)],beta,gam,k);

hts=fftshift(ifft(psihat));
htsre=real(hts);
htsim=imag(hts);
envel(k+1,:)=htsre.^2+htsim.^2;
end;
vetja=zeros(1,n);
for k=1:nk
vetja=vetja+lambda(k).*envel(k,:);
end;
vetja2=sqrt(vetja);
vector=sum(vetja2)^2/sum(vetja)/A(j)/2;