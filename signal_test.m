function [a1,a2,if1,if2,s1,s2,s] = signal_test(t)

a1 = exp(2*(1-t).^3 + 1.5*t.^4);
a2 = 1+ 5*t.^3 + 7*(1-t).^6;

phi1 = 50*t+30*t.^3-20*(1-t).^4;
phi2 = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));

if1 = 50+90*t.^2+80*(1-t).^3; 
if2 = 340+4*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))-28*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2)); 

s1 = a1.*exp(2*pi*1i*(phi1));
s2 = a2.*exp(2*pi*1i*(phi2));

s = s1+s2;
end
