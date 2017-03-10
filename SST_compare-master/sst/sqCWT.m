function [tfr, tfrsq, tfrtic, tfrsqtic] = sqCWT(t, x, freqlow, freqhigh, alpha, opts);

%
% Synchrosqueezing transform ver 0.5 (2015-03-09)
% You can find more information in 
%	http://sites.google.com/site/hautiengwu/
%
% Example: [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
%	time: 	time of the signal
%	xm: 	the signal to be analyzed
%	[lowfreq, highfreq]: the frequency range in the output time-frequency representation. For the sake of computational efficiency.
%	alpha:	the frequency resolution in the output time-frequency representation
%	opts:	parameters for the CWT analysis. See below
%	tfr/tfrtic:	the CWT and its scale tic
%	tfrsq/tfrsqtic: the SST-CWT and its frequency tic
%
% by Hau-tieng Wu v0.1 2011-06-20 (hauwu@math.princeton.edu)
%		  v0.2 2011-09-10
%		  v0.3 2012-03-03
%		  v0.4 2012-12-12
%		  v0.5 2015-03-09

	%% you can play with these 4 parameters, but the results might not be
	%% that different if the change is not crazy
Gamma = 1e-8;
nvoice = 32;
scale = 2;
oct = 1;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
noctave = floor(log2(length(x))) - oct;
dt = t(2) - t(1);

    
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %% Continuous wavelet transform 

[tfr, tfrtic] = CWT(t, x, oct, scale, nvoice, opts) ;
Dtfr = (-i/2/pi/dt)*[tfr(2:end,:) - tfr(1:end-1,:); tfr(end,:)-tfr(end-1,:)] ;

Dtfr((abs(tfr) < Gamma)) = NaN;
omega = Dtfr./tfr;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	%% Synchro-squeezing transform

[tfrsq, tfrsqtic] = SQ(tfr, omega, alpha, scale, nvoice, freqlow, freqhigh);
tfr = transpose(tfr) ;
tfrsq = transpose(tfrsq) ;    

end










%=====================================================================
	%% function for CWT-based SST
function [tfrsq, freq] = SQ(tfd, omega, alpha, scale, nvoice, freqlow, freqhigh);

omega = abs(omega);
[n, nscale] = size(tfd);

nalpha = floor((freqhigh - freqlow)./alpha);
tfrsq = zeros(n, nalpha);
freq = ([1:1:nalpha])*alpha + freqlow;
nfreq = length(freq);
	
for b = 1:n             %% Synchro-
    for kscale = 1: nscale       %% -Squeezing

        qscale = scale .* (2^(kscale/nvoice));

        if (isfinite(omega(b, kscale)) && (omega(b, kscale)>0))
            k = floor( ( omega(b,kscale) - freqlow )./ alpha )+1;

            if (isfinite(k) && (k > 0) && (k < nfreq-1))
	        ha = freq(k+1)-freq(k);
                tfrsq(b,k) = tfrsq(b,k) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
            end
        end
    end
end   
end
