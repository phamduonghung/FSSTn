% test_hm.m
%
% P. Flandrin & J. Xiao, 2005
%
% computes and displays Hermite multitaper (reassigned) spectrograms
% in the case of a synthetic sinusoidal FM signal embedded in white
% Gaussian noise
%
% calls  - tfrrsp_hm.m
%        - imageTF.m

% signal 

SNR = 10 ; % SNR
Nx = 256 ; % signal length
t = 1:Nx ;
[signal,iflaw] = fmsin(Nx,0.15,0.45,2*Nx,Nx/4,0.3,-1) ; 
x = sigmerge(signal,randn(size(signal)),SNR) ;

% analysis parameters

Nfft = Nx ; % number of FFT bins
M = 6 ; % number of Hermite tapers
Nh = 95 ; % length of Hermite functions

% computation

[S,RS,hat,tt] = tfrrsp_hm(x,t,Nfft,Nh,M,6) ;

SS = [] ; RSS = [] ;

for opt = 1:5 % loop on means (cf. mean_hm.m)
    
    for k = 1:M % loop on tapers
        
        [opt k]

        [Sm,RSm] = mean_hm(S(:,:,1:k),RS(:,:,1:k),hat(:,:,1:k),opt) ;
        SS(:,:,opt,k) = Sm ;
        RSS(:,:,opt,k) = RSm ;

    end

end

% display

col = 1 ; % colormap (cf. imageTF.m)

figure(1) % multitaper spectrograms

for opt = 1:4 % loop on means
    
    for k = 1:M % loop on # tapers
        
        subplot(5,M,k+M*(opt-1))
        imageTF(SS(1:Nfft/2,:,opt,k),40,col,0)
        xlabel('')
        ylabel('')

    end

end

subplot(5,M,1)
ylabel('arithm.')
subplot(5,M,M+1)
ylabel('geom.')
subplot(5,M,2*M+1)
ylabel('min.')
subplot(5,M,3*M+1)
ylabel('med.')

subplot(5,M,1)
title('1 taper')

for k = 2:M
    subplot(5,M,k)
    title([int2str(k),' tapers'])
end

figure(2) % multitaper reassigned spectrograms

for opt = 1:5 % loop on means
    
    for k = 1:M % loop on # tapers
        
        subplot(5,M,k+M*(opt-1))
        imageTF(RSS(1:Nfft/2,:,opt,k),40,col,0)
        xlabel('')
        ylabel('')

    end

end

subplot(5,M,1)
ylabel('arithm.')
subplot(5,M,M+1)
ylabel('geom.')
subplot(5,M,2*M+1)
ylabel('min.')
subplot(5,M,3*M+1)
ylabel('med.')
subplot(5,M,4*M+1)
ylabel('mean field')

subplot(5,M,1)
title('1 taper')

for k = 2:M
    subplot(5,M,k)
    title([int2str(k),' tapers'])
end