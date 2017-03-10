function colridge_hanford(h,C,me)

im = findall(h,'type','image');
ax = findall(h,'type','axes');
hold(ax);
cdat = get(im,'CData');

[nr N] = size(C);

mycol = {'red','green','magenta','blue','yellow',...
    'blue','green','red','cyan','magenta','yellow',...
    'blue','green','red','cyan','magenta','yellow'};

t = get(im,'XData');
f = get(im,'YData');
val = f;
for k=1:nr
    if me==2
        plot(t,f(C(k,:)),'color','r','Parent',ax,'LineWidth',1.5);
    else 
        plot(t,f(C(k,:)),'color','g','Parent',ax,'LineWidth',1.5);
    end    
end

end


