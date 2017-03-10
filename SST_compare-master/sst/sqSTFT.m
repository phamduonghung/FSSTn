function [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, lowFreq, highFreq, alpha, tDS, h, Dh);
%
% Synchrosqueezing modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
%
%	computes the STFT and its SST
%
%   Example:
%
%	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
%	x=cos(2*pi*(4*t+cos(t/2))) ;
%	[h, Dh] = hermf(71, 1, 6) ;
%		%% get the TF representation of the signal x with the frequency range
%		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
%	[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0.1, 0.4, 0.001, 1, h', Dh');
%	imageRTF(t, tfrsqtic*Hz, abs(tfrsq)) ;
%
%		%% the first version can be recovered by
%   [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0, 0.5, 0.5/length(x), 1, h', Dh');
%
%
%======================================
%	X     : analysed signal.
%	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
%	alpha : the resolution in the frequency axis
%	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
%	H     : frequency smoothing window, H(0) being forced to 1
%   DH    : differentiation of H	
%	TFR   : STFT
%	TFRSQ : synchrosqueezed STFT 
%
%	F. Auger, May-July 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

[xrow,xcol] = size(x) ;
t = [1:length(x)] ;
tLen = length(t(1:tDS:length(x))) ;
tfrsqtic = [lowFreq:alpha:highFreq] ;

N = length([lowFreq:alpha:highFreq]) ;
Lidx = ceil( (N/2)*(lowFreq/0.5) ) ; 
Hidx = floor( (N/2)*(highFreq/0.5) ) ; 
fLen = Hidx - Lidx + 1 ;
tfrsqtic = linspace(lowFreq, highFreq, fLen)' ;



%====================================================================
	%% check input signals
if (xcol~=1),
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) | (rem(tDS,1)) 
    error('tDS must be an integer value >= 1');
end; 

[hrow,hcol] = size(h); Lh = (hrow-1)/2; 
if (hcol~=1)|(rem(hrow,2)==0),
    error('H must be a smoothing window with odd length');
end;


%====================================================================
	%% run STFT and reassignment rule
tfr = zeros(N, tLen); 	% for h
tf3 = zeros(N, tLen);	% for Dh
tfrtic = linspace(0, 0.5, N/2)' ;

for tidx = 1:tLen,
    ti = t((tidx-1)*tDS+1); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    norm_h=norm(h(Lh+1+tau));
    tfr(indices,tidx) = x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
    tf3(indices,tidx) = x(ti+tau).*conj(Dh(Lh+1+tau)) /norm_h;
end ;

	%% FFT on columns
tfr = fft(tfr);
tf3 = fft(tf3);
avoid_warn = find(tfr~=0);
tf3(avoid_warn) = round(imag(N*tf3(avoid_warn)./tfr(avoid_warn)/(2.0*pi)));

tfr = tfr(1:N/2, :) ;
tfr3 = tfr(1:N/2, :) ;





%====================================================================
	%% run synchrosqueezing

tfrsq = zeros(fLen, tLen); 
Ex = mean(abs(x(min(t):max(t))).^2); 
Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex

for icol = 1:tLen,
    for jcol = 1: N/2,
  	if abs(tfr(jcol,icol)) > Threshold,

   	    jcolhat = jcol - tf3(jcol,icol);
   	    jcolhat = rem(rem(jcolhat-1,N)+N,N)+1;

   	    if (jcolhat < Hidx + 1) & (jcolhat >= Lidx)
		    tfrsq(jcolhat-Lidx+1,icol) = tfrsq(jcolhat-Lidx+1, icol) + tfr(jcol,icol) ;
	    end

  	end;
    end;
end;

