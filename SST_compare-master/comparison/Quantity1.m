% This code draws the synchrosqueezed transforms in Figure 9 (bottom) of
% H. Yang, Robustness Analysis of Synchrosqueezed Transforms, preprint,
% 2014.

clear all;
close all;
numTest = 20;

recMat = [];
NMvec = 0:0.2:4;
length(NMvec)
for cntNM = 1:length(NMvec)
    NM = NMvec(cntNM)
    recVec = zeros(1,7);
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
            t_sc = 1/2 + 2/8;
            lowfreq = 0; highfreq = 130;
            
            % apply the SSWPT
            T_f = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,highfreq,lowfreq,rad,is_cos,t_sc,red,epsl,res,0);
            EMDrec1 = EMDMat(real(T_f(end/2+1:end,:)),InstFreq,highfreq);
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % multitaper time-frequency reassignment
            
            % analysis parameters
            Nfft = Nx*5 ; % number of FFT bins
            M = 5; % number of Hermite tapers
            for Nh = 505
                [S,RS,hat,tt] = tfrrsp_hm(x,t,Nfft,Nh,M,6) ;
                SS = [] ; RSS = [] ;
                for opt = 1:2 % loop on means (cf. mean_hm.m)
                    for k = M % loop on tapers
                        [opt k]
                        [Sm,RSm] = mean_hm(S(:,:,1:k),RS(:,:,1:k),hat(:,:,1:k),opt) ;
                        SS(:,:,opt,k) = Sm ;
                        RSS(:,:,opt,k) = RSm ;
                    end
                end
                EMDrec2 = EMDMat(RSS(1:Nfft/8,:,1,M),InstFreq,highfreq);
                EMDrec3 = EMDMat(RSS(1:Nfft/8,:,2,M),InstFreq,highfreq);
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % multitaper time-frequency reassignment
            
            % analysis parameters
            Nfft = Nx*5 ; % number of FFT bins
            M = 10; % number of Hermite tapers
            for Nh = 505
                [S,RS,hat,tt] = tfrrsp_hm(x,t,Nfft,Nh,M,6) ;
                SS = [] ; RSS = [] ;
                for opt = 1:2 % loop on means (cf. mean_hm.m)
                    for k = M % loop on tapers
                        [opt k]
                        [Sm,RSm] = mean_hm(S(:,:,1:k),RS(:,:,1:k),hat(:,:,1:k),opt) ;
                        SS(:,:,opt,k) = Sm ;
                        RSS(:,:,opt,k) = RSm ;
                    end
                end
                EMDrec4 = EMDMat(RSS(1:Nfft/8,:,1,M),InstFreq,highfreq);
                EMDrec5 = EMDMat(RSS(1:Nfft/8,:,2,M),InstFreq,highfreq);
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % Synchrosqueezed wavelet transform (SSWT)
            
            % set parameters
            opts = struct();
            opts.motherwavelet = 'Cinfc' ;
            opts.CENTER = 1 ;
            opts.FWHM = 0.3 ;
            
            lowfreq = 0 ;
            highfreq = 130 ;
            Hz = 2^10 ;
            L = 1 ;
            time = [1/Hz:1/Hz:L]' ;
            N = length(time) ;
            alpha = 0.2 ;
            
            scrsz = get(0,'ScreenSize');
            xm = x;
            
            % apply the SSWT
            [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
            
            EMDrec6 = EMDMat(abs(tfrsq),InstFreq,highfreq);
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % Synchrosqueezed short-time Fourier transform (SSSTFT)
            
            % set parameters
            L = 1 ;
            N = length(x) ;
            Hz = N ;
            highFq = 130;
            opts.WinLen = 401;
            alpha = 0.025/Hz ;
            initstate(1) ;
            xm = x ;
            time = [1/Hz:1/Hz:L]' ;
            [h,Dh,t] = hermf(opts.WinLen, 1, 6) ;
            
            % apply the SSSTFT
            [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(xm, 0, 0.5*highFq/Hz*2, alpha, 1, h', Dh');
            
            EMDrec7 = EMDMat(abs(tfrsq),InstFreq,highfreq);
        end
        recVec = recVec + [EMDrec1 EMDrec2 EMDrec3 EMDrec4 EMDrec5 EMDrec6 EMDrec7]/numTest;
    end
    recMat = [recMat; recVec];
    head = sprintf('results/recQuan1_NTest_%d.mat',numTest);
    save(head,'recMat');
end

%%
if 1
    %head = sprintf('results/recQuan1_NTest_%d.mat',numTest);
    %load(head);
    [sm nplot] = size(recMat);
    
    figure;
    hold on;
    h = zeros(1, 5);
    x = 0:0.2:4;
    x = x(:);
    y = recMat(:,1);
    h(1) = plot(x,y,'-*k');
    y = recMat(:,4);
    h(3) = plot(x,y,'-^b');
    y = recMat(:,5);
    h(5) = plot(x,y,'-om');
    y = recMat(:,2);
    h(2) = plot(x,y,'-^r');
    y = recMat(:,3);
    h(4) = plot(x,y,'-oc');
    
    
    legend(h,'SSWPT10','MRSA5','MRSA10','MRSG5','MRSG10','Location','southeast');
    axis square;
    hold off;
    xlabel('\sigma^2');ylabel('EMD');
    
%     set(gca, 'FontSize', 18);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
%     str = sprintf('results/compRS%d',numTest);
%     print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

    figure;
    hold on;
    h = zeros(1, 3);
    x = 0:0.2:4;
    x = x(:);
    y = recMat(:,1);
    h(1) = plot(x,y,'-*k');
    y = recMat(:,6);
    h(2) = plot(x,y,'-^b');
    y = recMat(:,7);
    h(3) = plot(x,y,'-om');
    
    
    legend(h,'HSST','SSWT','SSSTFT','Location','southeast');
    axis square;
    hold off;
    xlabel('\sigma^2');ylabel('EMD');
    
%     set(gca, 'FontSize', 18);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
%     str = sprintf('results/compOTorg%d',numTest);
%     print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

% head = sprintf('rslt%d.mat',numTest);
% save(head);
