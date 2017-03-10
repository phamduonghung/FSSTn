function vector=nyttcoh444(x,x2,x3,A,beta,gam,nk,deltat,tnought)
% This function calculates the polarization using complex wavelets
% and also 
%
% x is the 1st data vector
% x2,x3 are the 2nd and 3rd data vector
% A is the scale vector 
% beta is the first morse wavelet parameter
% gam the second
% nk is the no of complex wavelets used
% deltat refers to the sample interval
% 
% It uses a convolution to perform the calculations: be cautious with small
% values of a, as these may give not unsubstantial leakage.
% Uses Park et al's expression for the covariance of the incremental vector
% ONLY graphs angular estimate and its uncertainty
% corrected sept 2001 for covariance measure

% stating the properties of the wavelets
spread=15;
amin=2*0.45*deltat;
lll=analyt6(beta,gam,24*pi,nk);

orient tall
bet=beta;
n=length(x);

%Setting up the necessary matrices
picture=zeros(length(A),n,nk);
picture2=zeros(length(A),n,nk);
picture21=zeros(length(A),n,nk);
picture3=zeros(length(A),n,3);
picture4=zeros(length(A),n,3);
picture5=zeros(length(A),n,3);
picture6=zeros(length(A),n);
picture7=zeros(length(A),n);

%calculating the wavelet transform
for k=0:nk-1
%picture(:,:,k+1)=wscal55b(x,A,beta,gam,k,deltat);
%picture2(:,:,k+1)=wscal55b(x2,A,beta,gam,k,deltat);
%picture21(:,:,k+1)=wscal55b(x3,A,beta,gam,k,deltat);
picture(:,:,k+1)=wscal55b(x,A,beta,gam,k,deltat);
picture2(:,:,k+1)=wscal55b(x2,A,beta,gam,k,deltat);
picture21(:,:,k+1)=wscal55b(x3,A,beta,gam,k,deltat);
end;

%calculating the singular value decomposition
for t=1:n
for a=1:length(A)

%setting up the M matrix (or its transpose)
for k=1:nk
M(1,k)=conj(picture(a,t,k));
M(2,k)=conj(picture2(a,t,k));
M(3,k)=conj(picture21(a,t,k));
end;

% Performing the svd
[U,ve,V]=svd(transpose(M),0);

%entrying the polarization singular values
picture3(a,t,1)=ve(1,1);
if length(ve)<2
picture3(a,t,2)=ve(1,1);
else
picture3(a,t,2)=ve(2,2);
end;
if length(ve)<3
picture3(a,t,3)=picture3(a,t,2);
else 
picture3(a,t,3)=ve(3,3);
end;


% Calculating the imaginary and real part of the first eigenvector, finding
% the proportion of the real magnitude square
sf=real(V(:,1));
sf2=imag(V(:,1));
vvv1=transpose(sf)*sf;
vvv2=transpose(sf2)*sf2;
vvv3(a,t)=vvv1/(vvv2+vvv1);
picture4(a,t,:)=sf;
picture5(a,t,:)=sf2;
picture3(a,t,:)=picture3(a,t,:).^2/(picture3(a,t,1)^2+picture3(a,t,2)^2+picture3(a,t,3)^2);
% Calculating angle and its variance
picture6(a,t)=atan2(imag(V(1,1)*conj(V(3,1))),real(V(1,1)*conj(V(3,1))));
ggj=(picture3(a,t,2)+picture3(a,t,3))/(2*(nk-1)*abs(V(1,1)*V(3,1))^2*picture3(a,t,1));
ggjs=V(:,2)*ctranspose(V(:,2))+V(:,3)*ctranspose(V(:,3));
ggjs2=ggjs(3,1);
ggjs1=ggjs(1,1);
ggjs3=ggjs(3,3);
%picture7(a,t)=ggj*(ggjs3*abs(V(1,1))^2+ggjs1*abs(V(3,1))^2-2*real(V(1,1)*imag(V(3,1))*ggjs2));
picture7(a,t)=0.5.*ggj*(ggjs3*abs(V(1,1))^2+ggjs1*abs(V(3,1))^2-2*real(V(1,1)*conj(V(3,1))*ggjs2));
picture7(a,t)=sqrt(picture7(a,t));
%picture8(a,t)=acos(dot(sf,sf2)/(norm(sf)*norm(sf2)));


end;
end;

%colormap(1-gray(256))
figure

imagesc(((tnought:n+tnought-1)).*deltat,(A*deltat),(picture3(:,:,1)))
colorbar('vert')
xlabel('time')
ylabel('frequency')
title('coherence')


clear picture5;
clear picture4;
clear picture3;


% Plotting the results



figure
subplot(2,1,1)
%colormap(1-gray(256))
%[C,h]=contourf(((tnought:n+tnought-1)).*deltat,A*deltat,abs(picture6),vcs2);
%[C,h]=contourf(((tnought:n+tnought-1)).*deltat,A*deltat,picture6,vcs2);

imagesc(((tnought:n+tnought-1)).*deltat,((A*deltat)),(picture6))
axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
colorbar('vert')
%caxis([0 3.2]);
%caxis([-3.3 3.3]);
%clabel(C,h)
xlabel('time');
ylabel('scale');
title('Phase Estimate')
% colorbar('vert')
%hold on
%plot(thm2,spread.*(thm2-tnought),'k-',thm3,spread.*(2*thm-thm3+tnought),'k-',[thm2 thm3],ones(1,8).*amin,'k-')
%axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
%hold off
%vcs2=[-3.3:0.3:3.3];
%vcs2=[0:0.15:2];
%subplot(3,1,2)
%colormap(1-gray(256))
%[C,h]=contourf(((tnought:n+tnought-1)).*deltat,A*deltat,abs(abs(picture6)-pi/2),vcs2);
%[C,h]=contourf(((tnought:n+tnought-1)).*deltat,A*deltat,picture6,vcs2);
%axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
%caxis([0 3.3]);
%caxis([-3.3 3.3]);
%clabel(C,h)
%xlabel('time');
%ylabel('scale');
%title('Absolute angle deviation')
%% colorbar('vert')
%hold on
%plot(thm2,spread.*(thm2-tnought),'k-',thm3,spread.*(2*thm-thm3+tnought),'k-',[thm2 thm3],ones(1,8).*amin,'k-')
%axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
%hold off
subplot(2,1,2)

imagesc(((tnought:n+tnought-1)).*deltat,(A*deltat),(picture7))
%imagesc(((tnought:n+tnought-1)).*deltat,A*deltat,picture7)
axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
caxis([0 1])
%clabel(C,h);
xlabel('time');
ylabel('scale');
title('Phase uncertainty')
colorbar('vert')
%hold on
%plot(thm2,spread.*(thm2-tnought),'k-',thm3,spread.*(2*thm-thm3+tnought),'k-',[thm2 thm3],ones(1,8).*amin,'k-')
%axis([tnought*deltat (n-1+tnought)*deltat A(1)*deltat A(length(A))*deltat]);
%hold off


