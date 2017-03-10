function [xpad, nup, n1, n2] = padsignal(x, padtype)
%
% written by Eugene Brevdo
%

[nup, n1, n2] = p2up(length(x));

xl = padarray(x(:), n1, padtype, 'pre');
xr = padarray(x(:), n2, padtype, 'post');

xpad = [xl(1:n1); x(:); xr(end-n2+1:end)];

end

function [up,n1,n2] = p2up(n)
    %up = 2^(1+round(log2(n+eps)));
    up = 2^(ceil(log2(n+eps)));
    n1 = floor((up-n)/2); n2 = n1;
    if (mod(2*n1+n,2)==1), n2 = n1 + 1; end
end
