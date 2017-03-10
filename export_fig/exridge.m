function [c,e] = exridge(Tx,lambda,beta)
% exridge : extracts  the ridge curve by maximising some energy. 
%   The result is an approximation computd by a greedy algorithm.
%   The algorithm uses several random initializations and then
%   forward/backward search.
%
% INPUTS:   
%   Tx : synchrosqueezed transform
%   lambda : the paremeter associated with first derivative penalization
%   beta   : the parameter associated with the second derivative penalization 
% OUTPUTS:  
%   c : vector containing the indexes of the ridge. length :size(Tx,2).
%   e : energy of the returned ridge, scalar.


%Et = log(abs(Tx)+eps^0.25);
Et     = abs(Tx).^2;
[na,N] = size(Tx);

% Parameters of the algorithm.
da = 10; % Maximal jump allowed
ng = 60; % Number of random initializations.

ccur = zeros(1,N);
c = ccur;    
e = -Inf;


for k = floor(linspace(N/(ng+1),N-N/(ng+1),ng))
    
    [ecur,idx] = max(Et(:,k-1));        
    ccur(k-1) = idx;
    [ecur,idx] = max(Et(:,k));        
    ccur(k) = idx;
  
    % forward step
    for b=k+1:N
     etmp = -Inf;
     for a=max(1,idx-da):min(na,idx+da)
      if Et(a,b)-lambda*(a-ccur(b-1))^2-beta*(a-2*ccur(b-1)+ccur(b-2))^2 > etmp
       etmp = Et(a,b)-lambda*(a-ccur(b-1))^2-beta*(a-2*ccur(b-1)+ccur(b-2))^2;
       idx = a;
      end
     end
     ccur(b)=idx;
     ecur = ecur + etmp;
    end

    % backward step
    idx = ccur(k);
    for b=k-1:-1:1
        etmp = -Inf;
        for a=max(1,idx-da):min(na,idx+da)
            if Et(a,b)-lambda*(a-ccur(b+1))^2-beta*(a-2*ccur(b+1)+ccur(b+2))^2 > etmp
                etmp = Et(a,b)-lambda*(a-ccur(b+1))^2-beta*(a-2*ccur(b+1)+ccur(b+2))^2;
                idx = a;
            end
        end
        ccur(b)=idx;
        ecur = ecur + etmp;
    end
    
    if ecur> e
        e = ecur;
        c = ccur;
    end
end

end
