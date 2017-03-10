function [MultiTaper, MultiTaperAll tfrsqtic] = sqSTFTmultitaper(x, lowFreq, highFreq, alpha, WinLen, dim, supp, Smooth) ;

%
% Usage: 
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFTmultitaper(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp)
%
% MT = 1: ordinary SST; MT > 1: ConceFT
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfrsq, MultiTaper, tfrsqtic] = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6);


N = length(x) ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	%% generate the window for short time Fourier transform (STFT)
[h, Dh, ~] = hermf(WinLen, dim, supp) ;


%=======================================

[~, ~, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, 1, h(1,:)', Dh(1,:)', Smooth, 0);

MultiTaperAll = zeros(size(tfrsq,1), size(tfrsq,2), dim) ;
MultiTaperAll(:,:,1) = tfrsq ;

fprintf(['MultiTaper STFT total (Smooth = ',num2str(Smooth),', Hemi = ',num2str(0),') : ',num2str(dim),'; now:     ']) ;
for ii = 2: dim

	fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;
	[~, ~, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, 1, h(ii,:)', Dh(ii,:)', Smooth, 0);

 	MultiTaperAll(:, :, ii) = tfrsq ;
end

MultiTaper = mean(MultiTaperAll, 3) ;
fprintf('\n') ;

end
