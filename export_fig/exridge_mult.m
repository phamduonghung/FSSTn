function [Cs, Es] = exridge_mult(Tx, nr, lambda, beta,clwin)
% Function exridge_mult : extracts the ridges of a multicomponent signal
% Inputs:
%   Tx : SQ transform of s
%   nr : number of ridges
%   lambda, beta : regularization parameters
%   clwin : frequency clearing window
% Outputs:
%   Cs Cs(:,j) is the j-th ridge location

if nargin<3
    lambda = 0.2;
    beta   = 0.4; 
    clwin  = 1;
elseif nargin<4
    clwin = 1;
end

Txs = Tx;

[na,N] = size(Txs);

Cs = zeros(N, nr);
Es = zeros(nr, 1);

for j=1:nr
    [Cs(:,j), Es(j)] = exridge(Txs,lambda,beta);

    % Remove this curve from the representation
    for b=1:N
        Txs(max(1,Cs(b,j)-clwin):min(na,Cs(b,j)+clwin),b)=0;
    end
    
end

Cs = Cs';

% sort the different ridges from HF to BF
tmp = sum(Cs,2);
[B idx] = sort(tmp,'descend');
Cs = Cs(idx,:);

end
