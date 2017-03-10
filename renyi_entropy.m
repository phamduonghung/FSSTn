function Ren=renyi_entropy(TFR,t,f,alpha)
%	Ren=renyi_entropy(TFR,t,f,alpha) calculates Renyi entropy from 2-D TFR
%   
%   Inputs:
%	TFR : (M,N) 2-D TFR function.
%	T : a time vector   (default : (1:N)).	
%	F : a frequency vector		(default : (1:M)).	
%	ALPHA : Renyi measure order	(default : 3).
%   
%   Outputs:
%   Ren=1/(1-ALPHA)*log2[Sum[TFR(Fi,Ti)^ALPHA dFi.dTi]]
%            Fi,Ti : Alpha-order Renyi entropy
%   ALPHA = 1: limit case, the outcomes will be Shannon entropy
%	Sha = - Sum[TFR(Fi,Ti)log2[TFR(Fi,Ti)]dFi.dTi]
%          Fi,Ti

if (nargin == 0),
 error('At least one parameter required');
end;

[M,N] = size(TFR);
if (nargin == 1),
 t=1:N; f=(1:M)'; alpha=3;
elseif (nargin == 2),
 f=(1:M)'; alpha=3;
elseif (nargin == 3),
 alpha=3;
end;

f=sort(f); %sort frequency vector in ascending order such that the first 
%row TFR must correspond to the lower frequencies

TFR = TFR./trapz(f,trapz(t,TFR,2)); 
% Normalisation TFR;
%trapz function is used to calculate 2D integral of %matrix TFR according
%to abscissa X and ordinate Y

if alpha == 1 % limit case case: Shannon entropy
 if (min(min(TFR))<0)
     error('distribution with negative values => alpha=1 not allowed');
 else
     Ren=-trapz(f,trapz(t,TFR.*log2(TFR+eps),2));
 end
else % Renyi entropy
    Ren=1/(1-alpha)*log2(trapz(f,trapz(t,TFR.^alpha,2))+eps);
end

