% hermf.m
%
% P. Flandrin & J. Xiao, 2005
%
% computes a set of orthonormal Hermite functions 
% (for use with tfrrsp_h.m) 
%
% input  : - N : number of points (must be odd) 
%          - M : maximum order 
%          - tm : half time support (>= 6 recommended) 
%
% output : - h : Hermite functions (MxN) 
%          - Dh : H' (MxN) 
%          - tt : time vector (1xN) 

function [h,Dh,tt] = hermf(N,M,tm) ; 

dt = 2*tm/(N-1) ; 
tt = linspace(-tm,tm,N) ; 
P = [] ; Htemp = [] ; DH = [] ; 
g = exp(-tt.^2/2) ; 

P(1,:) = ones(1,N) ; 
P(2,:) = 2*tt ; 

for k = 3 : M+1 
    
    P(k,:) = 2*tt.*P(k-1,:) - 2*(k-2)*P(k-2,:) ; 
    
end 

for k = 1:M+1 
    
    Htemp(k,:) = P(k,:).*g/sqrt(sqrt(pi)*2^(k-1)*gamma(k))*sqrt(dt) ; 
    
end 

h = Htemp(1:M,:) ; 

for k = 1:M 
    
    Dh(k,:) = (tt.*Htemp(k,:) - sqrt(2*(k))*Htemp(k+1,:))*dt ; 
    
end 
