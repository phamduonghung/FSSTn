% This code draws the synchrosqueezed transforms in Figure 9 (top-left) of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

clear all;
numTest = 20;

redNum = [1 5:5:50];
recMat = [];
for SNR = 40:-10:10
    recVec = zeros(1,length(redNum));
    for cntt = 1:numTest
        for cntred2 = 1:length(redNum)
            %% set up data
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
                
                is_real = 0;
                NM = 1/(10^(SNR/10));
                ns = NM*(randn(1,N)+randn(1,N));
                
                if is_real
                    ns = real(ns);
                    f2 = real(f2);
                end
                x = f2.' + ns.';
                signal = f2.';
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if 1
                % set up parameters
                res = 0.2;
                NG = N;
                is_unif = 1;
                typeNUFFT = 1;
                is_cos = 1;
                epsl = 1e-2;
                xo = x;
                rad = 1.3;
                fff = x.';
                red = redNum(cntred2);
                t_sc = 1/2 + 2/8;
                lowfreq = 0; highfreq = 130;
                
                % apply the SSWPT
                T_f = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,highfreq,lowfreq,rad,is_cos,t_sc,red,epsl,res,0);
                EMDrec = EMDMat(real(T_f(end/2+1:end,:)),InstFreq,highfreq);
            end
            recVec(cntred2) = recVec(cntred2) + EMDrec/numTest;
        end
    end
    recMat = [recMat; recVec];
    head = sprintf('results/recQuan2_NTest_%d.mat',numTest);
    save(head,'recMat');
end

if 1
    %head = sprintf('results/recQuan2_NTest_%d.mat',numTest);
    %load(head);
    [nplot sn] = size(recMat);
    
    pic = figure;
    hold on;
    h = zeros(1, nplot);
    x = [1 5:5:50];
    y = recMat(1,:);
    h(1) = plot(x,y,'-ok');
    y = (recMat(2,:));
    h(2) = plot(x,y,'-^k');
    y = (recMat(3,:));
    h(3) = plot(x,y,'-ob');
    y = (recMat(4,:));
    h(4) = plot(x,y,'-^b');
    
    
    legend(h,'SNR=40','SNR=30','SNR=20','SNR=10','Location','northeast');
    axis square;
    hold off;
    xlabel('Reduncancy');ylabel('EMD');
%     
%     set(gca, 'FontSize', 18);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
%     str = sprintf('results/compRed%d',numTest);
%     print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
%     head = sprintf('results/compRed%d.fig',numTest);
    %saveas(pic,head);
end
