function [STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn_test(s,gamma,sigma,ft,bt)
%% sstn : computes the STFT of a signal and different versions of synchrosqueezing/reassignment.
%   Uses a Gaussian window.
%
% INPUTS:  
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:   
%   STFT: the short-time Fourier transform
%   SST1: standard synchrosqueezing
%   SST2: vertical second-order synchrosqueezing
%   SST3: vertical third-order synchrosqueezing
%   SST4: vertical fourth-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   tau2: second-order phase delay (horizontal reassignment operator)
%   tau3: third-order phase delay (horizontal reassignment operator)
%   omega2: second-order instantaneous frequency
%   omega3: third-order instantaneous frequency


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
sz=zeros(n,1);
% Padding
sleft = flipud(conj(sz(2:n/2+1)));
sright = flipud(sz(end-n/2:end-1));
x = [sleft; s ; sright];
clear xleft xright;

% Window definition
t = -0.5:1/n:0.5-1/n;t=t';
g =  1/sigma*exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'
%gpp = (-2*pi/sigma^2+4*pi^2/sigma^4*t.^2) .* g; % g''

% Initialization
STFT = zeros(neta,nb);
SST1 = zeros(neta,nb);
SST2 = zeros(neta,nb);
SST3 = zeros(neta,nb);
SST4 = zeros(neta,nb);
omega = zeros(neta,nb);
tau2 = zeros(neta,nb);
tau3 = zeros(neta,nb);
tau4 = zeros(neta,nb);
omega2 = zeros(neta,nb);
omega3 = zeros(neta,nb);
omega4 = zeros(neta,nb);
phi22p = zeros(neta,nb);
%phi2p = zeros(neta,nb);
phi23p = zeros(neta,nb);
phi33p = zeros(neta,nb);
%phi3p = zeros(neta,nb);
phi24p = zeros(neta,nb);
phi34p = zeros(neta,nb);
phi44p = zeros(neta,nb);
% = zeros(neta,nb);
vg = zeros(neta,7);
vgp = zeros(neta,5);
Y = zeros(neta,4,4);

%% Computes STFT and reassignment operators
for b=1:nb
	
    % STFT, window t^n*g
    for i = 0:7
        tmp = (fft(x(bt(b):bt(b)+n-1).*(t.^i).*g))/n;
        vg(:,i+1) = tmp(ft);
    end 
       
      %% STFT, window gx^np
    for i = 0:5
        tmp = fft(x(bt(b):bt(b)+n-1).*(t.^i).*gp)/n;
        vgp(:,i+1) = tmp(ft);
    end
    
    %% second-order operator tau
    tau2(:,b) = vg(:,2)./vg(:,1);
    
    % third order operator tau
    tau3(:,b) = vg(:,3)./vg(:,1);
    
    % four order operator tau
    tau4(:,b) = vg(:,4)./vg(:,1);
        
    %% Y expressions
    for i = 1:7
        for j = 1:7
            if i>=j
                Y(:,i,j) = vg(:,1).*vg(:,i+1) - vg(:,j).*vg(:,i-j+2);
            end
        end
    end  
    %% W expressions
    W2 = 1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));
    W3 = 1/2/1i/pi*(2*vg(:,1).*vg(:,2)+vg(:,1).*vgp(:,3)-vg(:,3).*vgp(:,1));
    W4 = 1/2/1i/pi*(2*vg(:,1).*vg(:,3)+2*vg(:,2).^2+vg(:,1).*vgp(:,4) - vg(:,4).*vgp(:,1)+vg(:,2).*vgp(:,3) - vg(:,3).*vgp(:,2));
    
    %% operator omega
    omega(:,b) = (ft-1)'-real(vgp(:,1)/2/1i/pi./vg(:,1));

    %% operator hat p: estimations of frequency modulation  
    %SST2
    phi22p(:,b) = W2./Y(:,2,2);
    omega2(:,b) = omega(:,b) + real(phi22p(:,b).*tau2(:,b));
         
    %SST3
    phi33p(:,b) = (W3.*Y(:,2,2)-W2.*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));
    phi23p(:,b) = W2./Y(:,2,2) - phi33p(:,b).*Y(:,3,2)./Y(:,2,2);
    omega3(:,b) = omega(:,b) + real(phi23p(:,b).*tau2(:,b))+ real(phi33p(:,b).*tau3(:,b));
       
    %SST4    
    phi44p(:,b) =((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*W4-(W3.*Y(:,2,2)-W2.*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(W3.*Y(:,3,2)-W2.*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)))...
        ./((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*(Y(:,6,4)+Y(:,6,3)-Y(:,6,2))-(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(Y(:,5,3).*Y(:,3,2)-Y(:,4,2).*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)));
    
    phi34p(:,b) = (W3.*Y(:,2,2)-W2.*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3))-phi44p(:,b).*(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));
    phi24p(:,b) = W2./Y(:,2,2) - phi34p(:,b).*Y(:,3,2)./Y(:,2,2) - phi44p(:,b).*Y(:,4,2)./Y(:,2,2);
    omega4(:,b) = omega(:,b) + real(phi24p(:,b).*tau2(:,b))+ real(phi34p(:,b).*tau3(:,b))+ real(phi44p(:,b).*tau4(:,b));
     
    % Storing STFT
    STFT(:,b) = vg(:,1).* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
        
end

%% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
           
%%%%%%%%%%%%%%%%%%%%%%%%%%SST1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            k = 1+round(omega(eta,b));          
            if k>=1 && k<=neta               
             SST1(k,b)  = SST1(k,b) + STFT(eta,b);
            end
            
%%%%%%%%%%%%%%%%%%%%SST2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            k = 1+round(omega2(eta,b));
             if k>=1 && k<=neta
               SST2(k,b)  = SST2(k,b) + STFT(eta,b);
             end
             
%%%%%%%%%%%%%%%%%%%%%SST3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
             k = 1+floor(omega3(eta,b));
             if k>=1 && k<=neta
               SST3(k,b)  = SST3(k,b) + STFT(eta,b);
             end
             
%%%%%%%%%%%%%%%%%%%%%SST4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             k = 1+floor(omega4(eta,b));
             if k>=1 && k<=neta
               SST4(k,b)  = SST4(k,b) + STFT(eta,b);
             end

        end
    end
end