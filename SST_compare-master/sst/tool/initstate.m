function initstate(state)
%
% Initialize the state of the random numbers generator.
% This enables reproducing the generated "random" numbers in different
% calls to rand and randn.
%
% Yoel Shkolnisky, January 2007.

if nargin<1
    state = 3498;
end

rand('state', state);
randn('state', state);

