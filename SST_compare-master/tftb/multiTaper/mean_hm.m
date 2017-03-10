% mean_hm.m
%
% J. Xiao & P. Flandrin, 2005
%
% computes various means of multitaper (reassigned) spectrograms
% (for use with tfrrsp_hm.m)
%
% input  - S : spectrograms (Nfft x Nx x M)
%        - RS : spectrograms (Nfft x Nx x M)
%        - hat : reassignment vectors fields (Nfft x Nx x M)
%        - opt : option for type of means, 1 = arithmetic (default),
%          2 = geometric, 3 = min, 4 = median, 5 = mean spectrogram
%          reassigned with median reassignment vector field
%
% output - Sm : mean spectrogram (Nfft x Nx)
%        - RSm : mean reassigned spectrogram (Nfft x Nx)

function [Sm,RSm] = mean_hm(S,RS,hat,opt) ;

if opt == 1
    
    Sm = mean(S,3) ;
    RSm = mean(RS,3) ;
    
elseif opt == 2
    
    prec = 1e-12 ;
    Sm = exp(mean(log(max(S,prec)),3)) ;
    RSm = exp(mean(log(max(RS,prec)),3)) ;
    
elseif opt == 3
    
    Sm = min(S,[ ],3) ;
    RSm = min(RS,[ ],3) ;

elseif opt == 4
    
    Sm = median(S,3) ;
    RSm = median(RS,3) ;

elseif opt == 5
    
    [Nfft,N,C] = size(S) ;
    
    hatm = mean(hat,3) ;
    
	Sm = mean(S,3) ; 
    tfr = Sm ;
	tcol = N ;
	rtfr= zeros(Nfft,tcol) ; 
	Threshold = 1.0e-10*max(max(tfr)) ;
    
	for icol=1:tcol

        for jcol=1:Nfft
      
            if abs(tfr(jcol,icol))>Threshold,
      
                icolhat = round(imag(hatm(jcol,icol))) ;
                icolhat = min(max(icolhat,1),tcol) ;
                jcolhat = mod(round(real(hatm(jcol,icol))),Nfft) ;
                jcolhat = min(max(jcolhat,1),Nfft) ; 
                rtfr(jcolhat,icolhat) = rtfr(jcolhat,icolhat) + tfr(jcol,icol) ;
         
            else

                rtfr(jcol,icol) = rtfr(jcol,icol) + tfr(jcol,icol) ;
      
            end
     
        end
	
    end
    
    RSm = rtfr ;

end
