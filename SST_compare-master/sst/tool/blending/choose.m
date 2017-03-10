function k= choose(n, m)
%
% choose(n,m) is the number of ways of choosing m objects from n distinct
% objects.  The simplest definition of choose(n,m) is n! / (m! * (n-m)!).
% The algorithm used here is somewhat less susceptible to overflow, and is
% faster than Matlab's builtin nchoosek function.
%
% Dr. Phillip M. Feldman   21 April 2006


% Check input arguments:

if (nargin ~= 2)
   error('choose requires 2 arguments.');
end

if (m < 0 | m > n)
   error('m (second argument) must be between 0 and n, inclusive.');
end

% The simplest definition of choose(n,m) is n! / (m! * (n-m)!), but the
% following algorithm is somewhat less susceptible to overflow.

if (m >= n-m)
   k= prod(m+1:n) / prod(2:n-m);

else
   k= prod(n-m+1:n) / prod(2:m);
end
end
