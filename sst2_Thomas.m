function [STFT SST RSTFT OSST VSST GAR omega tau omega2] = sst2(s,gamma,sigma,ft,bt)
%% sst2 : computes the STFT of a signal and different versions of synchrosqueezing/reassignment.
%   Uses a Gaussian window.
%
% INPUTS:  
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
% INPUTS:   
%   STFT: the short-time Fourier transform
%   SST: standard synchrosqueezing
%   RSTFT: reassigned STFT
%   OSST: oblique synchrosqueezing
%   VSST: vertical second-order synchrosqueezing
%   GAR: modified complex reassigned STFT [2]
%   omega: instantaneous frequency (vertical reassignment operator)
%   tau: phase delay (horizontal reassignment operator)
%   omega2: second-order instantaneous frequency


% checking length of signal
n = length(s);
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);


% Optional parameters
if nargin<5
   ft = 1:n/2;
   bt = 1:n;
end
nb = length(bt);
neta = length(ft);


% Padding
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
clear xleft xright;

% Window definition
t = linspace(-0.5,0.5,n);t=t';
g = exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'
gpp = (-2*pi/sigma^2+4*pi^2/sigma^4*t.^2) .* g; % g''

% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
RSTFT = zeros(neta,nb);
OSST = zeros(neta,nb);
VSST = zeros(neta,nb);
GAR = zeros(neta,nb);
omega = zeros(neta,nb);
tau = zeros(neta,nb);
omega2 = zeros(neta,nb);
phipp = zeros(neta,nb);

%% Computes STFT and reassignment operators
df = ft(2)-ft(1);
db = bt(2)-bt(1);
for b=1:nb
	% STFT, window g
    tmp = (fft(x(bt(b):bt(b)+n-1).*g))/n;
    vg = tmp(ft);
    
    % STFT, window gp
    tmp = fft(x(bt(b):bt(b)+n-1).*gp)/n;
    vgp = tmp(ft);
    
    % operator omega
    omega(:,b) = (ft-1)'-real(tmp(ft)/2/1i/pi./vg);
    
    % STFT, window xg
    tmp = fft(x(bt(b):bt(b)+n-1).*t.*g)/n;
    vxg = n*tmp(ft);
    
    % operator tau
    tau(:,b) = n*real(tmp(ft)./vg);

    % STFT, window xxg
    tmp = fft(x(bt(b):bt(b)+n-1).*t.^2.*g)/n;
    vxxg = n^2*tmp(ft);
    
    % STFT, window gpp
    tmp = fft(x(bt(b):bt(b)+n-1).*gpp)/n;
    vgpp = tmp(ft);
    
    % STFT, window vxgp
    tmp = fft(x(bt(b):bt(b)+n-1).*t.*gp)/n;
    vxgp = n*tmp(ft);
        
    % operator hat q: estimation of frequency modulation
    %phipp(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg));
	phipp(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vg.^2+vxg.*vgp-vxgp.*vg));
    
    %phipp2(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./vg.^2)./real((vxg.*vgp-vxgp.*vg)./vg.^2);
    
    %phipp2(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./vg.^2)./real((vxg.*vgp-vxgp.*vg)./vg.^2);
    
    % Second-order instantaneous frequency
    omega2(:,b) = omega(:,b) - phipp(:,b).*tau(:,b);
    
    % Storing STFT
    STFT(:,b) = vg.* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
        
end

%imagesc(abs(phipp-phipp2));sum(sum(abs(phipp-phipp2)));error('ah')
%phipp(abs(STFT)<0.00001)=NaN;phipp2(abs(STFT)<0.00001)=NaN;figure();imagesc(phipp);figure();imagesc(phipp2);set(findall(0,'type','axes'),'clim',[-1 1]);error('bouh');

%% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
            k = 1+round((omega(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                % Vertical reassignment: SST 
                SST(k,b) = SST(k,b) + STFT(eta,b);
                %l = 1+round((tau(eta,b)-bt(1)+1)/db);
                l = round((tau(eta,b))/db);
                if b+l>=1 && b+l<=nb
                    % Oblique reassignment: RSTFT, OSST, GAR, 
                    RSTFT(k,b+l) = RSTFT(k,b+l) + abs(STFT(eta,b));
                    OSST(k,b+l) = OSST(k,b+l) + STFT(eta,b)*exp(1i*pi*(2*omega(eta,b)-phipp(k,b+l)*l*db)*l*db/n);
                    GAR(k,b+l) = GAR(k,b+l) + STFT(eta,b)*exp(1i*pi*(omega(eta,b)+(ft(eta)-1))*l*db/n);
                end
            end
            k = 1+round((omega2(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                % second-order Vertical reassignment: VSST
                VSST(k,b) = VSST(k,b) + STFT(eta,b);
            end
        end
    end
end

STFT = df*STFT;
SST = df*SST;
VSST = df*VSST;
OSST = df*OSST;
RSTFT = df*RSTFT;
GAR = df*GAR;
