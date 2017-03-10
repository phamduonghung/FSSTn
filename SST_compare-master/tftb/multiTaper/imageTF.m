% imageTF.m
%
% P. Flandrin, 2005
%
% displays a time-frequency distribution
%
% input - M : time-frequency matrix
%       - dyn : dynamic range (in dB)
%       - c : choice of colormap, 1 = jet (default), else = b&w
%       - b : colorbar, 1 = yes (default), else = no

function imageTF(M,dyn,c,b) ;

if (nargin == 2)
    
    c = 1 ;
    b = 1 ;
    
end

if (nargin == 3)
    
    b = 1 ;
    
end

M = M/max(max(M)) ;
M = max(M,10^(-dyn/10)) ;
imagesc(flipud(10*log10(M)),[-dyn 0])
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')

if c == 1
    
    colormap(jet)
    
else
    
    colormap(flipud(gray))
    
end

if b == 1
    
    colorbar
    
end
