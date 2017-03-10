function vector=wscal55b(x,A,beta,gam,k,deltat)
% This function calculates the cwt using the complex morsewavelets
% where
%
% x is the data vector
% A is the scale vector 
% beta is the first morse wavelet parameter
% gam the second
% k is the kth morsewavelets used
% deltat refers to the sample interval
% 
% It uses a convolution to perform the calculations

n=length(x);
q=size(x);
if q(1)>1
x=transpose(x);
end;
picture=zeros(length(A),n);
xhat=fft(x);
for j=1:length(A)
psihat(1:ceil(n/2))=sqrt(deltat.*abs(A(j))).*wwhat(A(j).*[(0:ceil(n/2)-1)/(n)],beta,gam,k);
% psihat(ceil(n/2)+1:n)=sqrt(deltat.*abs(A(j))).*wwhat(A(j).*[(ceil(n/2)+1:n)/(n)],beta,gam,k);
psihat(ceil(n/2)+1:n)=sqrt(deltat.*abs(A(j))).*wwhat(A(j).*[fliplr(-(1:floor(n/2)))/(n)],beta,gam,k);
what=xhat.*psihat;
vect2=ifft(what);
picture(j,:)=vect2;
end;
vector=picture;

