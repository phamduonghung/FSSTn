% This code draws the synchrosqueezed transforms in Figure 9 (top-right) of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

clear all;
close all;
numTest = 20;

recMat = [];
NMvec = 0:0.2:4;
t_scvec = 0.5:0.125:1;
length(NMvec)
for cntNM = 1:length(NMvec)
    NM = NMvec(cntNM)
    recVec = zeros(1,length(t_scvec));
    for cntt = 1:numTest
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
            ns = NM*(randn(1,N)+randn(1,N));
            
            if is_real
                ns = real(ns);
                f2 = real(f2);
            end
            x = f2.' + ns.';
            signal = f2.';
            SNR = 10*log10((var(f2(:)))/var(ns(:)))
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % SSWPT
            recEMD = []
            for cntt = 1:length(t_scvec)
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
                red = 10;
                t_sc = t_scvec(cntt);
                lowfreq = 0; highfreq = 130;
                
                % apply the SSWPT
                T_f = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,highfreq,lowfreq,rad,is_cos,t_sc,red,epsl,res,0);
                recEMD = [recEMD EMDMat(real(T_f(end/2+1:end,:)),InstFreq,highfreq)];
            end
        end
        
        recVec = recVec + recEMD/numTest;
    end
    recMat = [recMat; recVec];
    head = sprintf('results/recQuan3_NTest_%d.mat',numTest);
    save(head,'recMat');
end
%%
if 1
    %head = sprintf('results/recQuan3_NTest_%d.mat',numTest);
    %load(head);
    [sm nplot] = size(recMat);
    
    pic = figure;
    hold on;
    h = zeros(1, 5);
    x = NMvec;
    x = x(:);
    y = recMat(:,1);
    h(1) = plot(x,y,'-*k');
    y = recMat(:,2);
    h(2) = plot(x,y,'-^b');
    y = recMat(:,3);
    h(3) = plot(x,y,'-om');
    y = recMat(:,4);
    h(4) = plot(x,y,'-^r');
    y = recMat(:,5);
    h(5) = plot(x,y,'-oc');
    
    
    legend(h,'0.5','0.625','0.75','0.875','1','Location','southeast');
    axis square;
    hold off;
    xlabel('\sigma^2');ylabel('EMD');
    
    %     set(gca, 'FontSize', 18);
    %     b=get(gca);
    %     set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    %     str = sprintf('results/compt_sc%d',numTest);
    %     print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    %     head = sprintf('results/compt_sc%d.fig',numTest);
    %     saveas(pic,head);
end