function OTD = slicedOT(mu, nu, p)
%SLICEDOA Sliced Optimal-Transports between time-frequency images
%   compute OT on each fixed time slice, then take the average
%   mu, nu: input source/target measures
%   p: compute the p-wasserstein distance for each slice

if ~all(size(mu)==size(nu))
    error('dimensions do not match');
end

% NF = size(mu, 1);
NT = size(mu, 2);

if (nargin < 3)
    p = 1;
end

OTDs = zeros(NT,1);

if (p == 1)
    %%% In dimension 1, Wasserstein-1 distance is the same as L1 distance
    for j=1:NT
        s = mu(:,j); s = s/sum(s);
        t = nu(:,j); t = t/sum(t);        
        OTDs(j) = mean(abs(cumsum(s)-cumsum(t)));
    end
else
    lattice = 0:0.01:1;
    %%% In dimension 1, Wasserstein-p (p>1) distance is the L1 distance
    %%% between inverses of cumulative distribution functions
    for j=1:NT
        s = mu(:,j); s = s/sum(s);
        t = nu(:,j); t = t/sum(t);
        cdf_s = cumsum(s);
        cdf_t = cumsum(t);
%         cdf_s = interp1(1:length(cumsum(s)),cumsum(s),1:0.5:length(cumsum(s)));
%         cdf_t = interp1(1:length(cumsum(t)),cumsum(t),1:0.5:length(cumsum(t)));
        is = invCDF(cdf_s, lattice);
        it = invCDF(cdf_t, lattice);
        OTDs(j) = (mean(abs(is-it).^p))^(1/p);
    end
end

if sum(isnan(OTDs)); fprintf(['\t\t\t ****', num2str(sum(isnan(OTDs))),'\n']) ; end
OTD = mean(OTDs(find(~isnan(OTDs))));


end

function is = invCDF(s,l)
%%% s is a cdf defined equi-distantly on [1,2,\cdots,length(s)]
%%% l is strictly increasing and equi-distant points from 0 to 1
%%% l(1) = 0 and l(end) = 1
if (l(1) ~= 0)
    error('second argument should start with 0');
end
if (l(end) ~= 1)
    error('second argument should end with 1');
end
is = zeros(1,length(l)-1);
for j=1:(length(l)-1)
    is(j) = find(s>l(j), 1)/length(s);
end
end
