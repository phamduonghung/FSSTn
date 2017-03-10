% imageSQ.m
%
% Display time-frequency result of the Synchrosqueezing transform and others
% function imageSQ(t, ytic, M) ;
%
function imageSQ(t, ytic, M, Qv) ;
if nargin < 4, Qv = 1; end;

fz = 20;

S = size(M);
Q = M(:);
	
fprintf(['\n\t\t\t ** Smallest value = ',num2str(min(M(:))),'\n']) ;
fprintf(['\t\t\t ** Maximal value = ',num2str(max(M(:))),'\n']) ;
fprintf(['\t\t\t ** 99.8%% value = ',num2str(quantile(M(:),0.998)),'\n']) ;
fprintf(['\t\t\t ** 0.2%% value = ',num2str(quantile(M(:),0.002)),'\n']) ;

	% truncate the upper bound
q = quantile(Q, Qv);
M(find(M>q)) = q;
M = M ./ q ;

	% truncate the lower bound
m = quantile(Q, 0.002) ;
M(find(M<m)) = m ;


imagesc(t, ytic, M)
axis xy ;
set(gca, 'fontsize', fz);
