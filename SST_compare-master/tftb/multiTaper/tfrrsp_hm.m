% tfrrsp_hm.m
%
% J. Xiao & P. Flandrin, 2005
%
% computes Hermite multitaper (reassigned) spectrograms
%
% input  - x : signal
%        - t : time instants of analysis
%        - Nfft : number of frequency bins for FFT
%        - Nh : number of points for Hermite tapers (must be odd) 
%        - M : maximum order 
%        - tm : half time support (>= 6 recommended)
%
% output - S : spectrograms (Nfft x length(t) x M)
%        - RS : reassigned spectrograms (Nfft x length(t) x M)
%        - hat : reassigned vector fields (Nfft x length(t) x M)
%        - tt : Hermite support (1 x Nh)
%
% calls  - tfrrsp_h.m
%        - hermf.m

function [S,RS,hat,tt] = tfrrsp_hm(x,t,Nfft,Nh,M,tm) ;

[h,Dh,tt] = hermf(Nh,M,tm) ;

S = zeros(Nfft,length(t),M) ; 
RS = zeros(Nfft,length(t),M) ; 
hat = zeros(Nfft,length(t),M) ; 

for k = 1:M
    
    [spt,rspt,hatt] = tfrrsp_h(x,t,Nfft,h(k,:)',Dh(k,:)') ;
    S(:,:,k) = spt ;
    RS(:,:,k) = rspt ;
    hat(:,:,k) = hatt ;
    
end