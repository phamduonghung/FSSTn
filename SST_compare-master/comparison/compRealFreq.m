% This code draw the real instantaneous frequency in Figure 1 of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.


if(1)
    N = 2^10;
    Nx = N;
    t=1:Nx;
    fprintf('sampling rate is %d per second\n',N);
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.05;
    F2 = 30;
    yy = x +amp*cos(2*pi*x);
    f2 = exp(2*pi*i*F2*yy);
    InstFreq = F2*(1-amp*sin(2*pi*x)*2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sm = 652/4; sn = 1024;
temp = [];
for cnt = 1:sn
    distB = zeros(sm,1);
    distB(round(sm*InstFreq(cnt)/130)) = 1;
    temp = [temp distB];
end
figure;imagesc([0,1],[0,130],temp);axis square;         axis xy;
xlabel('Time') ; ylabel('Freq') ; %title('SST (display in linear scale)') ;
% set(gca, 'FontSize', 18);
% b=get(gca);
% set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
% str = 'results/nsReal';
% print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

