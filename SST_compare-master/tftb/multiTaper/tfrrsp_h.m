% tfrrsp_h.m
%
% P. Flandrin & J. Xiao, 2005, adapted from F. Auger's tfrrsp.m 
% in the Time-Frequency ToolBox (available at http://tftb.nongnu.org)
%
% computes a spectrogram and its reassigned version
% in the specific case of Hermite functions (for use
% with tfrrsp_hm.m and hermf.m)
%
% input  - x : signal
%        - t : time instants of analysis
%        - Nfft : number of frequency bins for FFT
%        - h : Hermite function (output of hermf.m) 
%        - Dh : derivative of h (output of hermf.m) 
%
% output - S : spectrogram (Nfft x length(t))
%        - RS : reassigned spectrogram (Nfft x length(t))
%        - hat : reassigned vector field (Nfft x length(t))

function [S,RS,hat] = tfrrsp_h(x,t,Nfft,h,Dh) ;

[xrow,xcol] = size(x) ;
hlength = floor(Nfft/4) ; hlength = hlength+1-rem(hlength,2) ;
[trow,tcol] = size(t) ; 
[hrow,hcol] = size(h) ; Lh=(hrow-1)/2 ; 

if (tcol==1)  
    
    Dt = 1 ;    
    
else    
    
    Deltat = t(2:tcol)-t(1:tcol-1) ; 
    Mini = min(Deltat) ; Maxi=max(Deltat);
    
    if (Mini~=Maxi)
        
        error('The time instants must be regularly sampled.') ;
        
    else
        
        Dt = Mini ;
        
    end 
    
    clear Deltat Mini Maxi
    
end

S = zeros(Nfft,tcol) ; tf2 = zeros(Nfft,tcol); tf3 = zeros (Nfft,tcol);

Th = h.*[-Lh:Lh]'; 

for icol = 1:tcol 
    
    ti = t(icol) ; 
    tau = -min([round(Nfft/2)-1,Lh,ti-1]):min([round(Nfft/2)-1,Lh,xrow-ti]) ;
    indices = rem(Nfft+tau,Nfft)+1 ;
    norm_h = norm(h(Lh+1+tau)) ;
    S(indices,icol) = x(ti+tau).*conj( h(Lh+1+tau))/norm_h ;
    tf2(indices,icol) = x(ti+tau).*conj(Th(Lh+1+tau))/norm_h ;
    tf3(indices,icol) = x(ti+tau).*conj(Dh(Lh+1+tau))/norm_h ;
    
end

S = fft(S) ; tf2 = fft(tf2) ; tf3=fft(tf3) ;
avoid_warn = find(S~=0) ;
tf2(avoid_warn) = round(real(tf2(avoid_warn)./S(avoid_warn)/Dt)) ;
tf3(avoid_warn) = round(imag(Nfft*tf3(avoid_warn)./S(avoid_warn)/(2*pi))); 
S = abs(S).^2 ;

RS = zeros(Nfft,tcol) ; 
Ex = mean(abs(x(min(t):max(t))).^2) ; Threshold=1.0e-6*Ex ;

for icol = 1:tcol
    
    for jcol = 1:Nfft
        
        if abs(S(jcol,icol))>Threshold
            
            icolhat = icol + tf2(jcol,icol) ;
            icolhat = min(max(icolhat,1),tcol) ;
            jcolhat = jcol - tf3(jcol,icol) ;
            jcolhat = rem(rem(jcolhat-1,Nfft)+Nfft,Nfft)+1 ;
            RS(jcolhat,icolhat) = RS(jcolhat,icolhat) + S(jcol,icol) ;
            tf2(jcol,icol) = jcolhat + j*icolhat; 
            
        else
            
            tf2(jcol,icol) = inf*(1+j) ;
            RS(jcol,icol) = RS(jcol,icol) + S(jcol,icol) ;
            
        end
        
    end
    
end

hat = tf2 ;