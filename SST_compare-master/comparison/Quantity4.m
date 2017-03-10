% This code compare the highly redundant synchrosqueezed wave packet
% transform with the ConceFT transform.
%
% By Haizhao Yang

clear all;
close all;
numTest = 20;
NMvec = 0:0.2:1;

% This code is very slow. You can plot the results precomputed or run this
% code to compute the results again.

isCompute = 0; % whether compute the results again, 0: yes, 1: no

if isCompute
recMat = [];
length(NMvec)
for cntNM = 1:length(NMvec)
    NM = NMvec(cntNM)
    recVec = zeros(1,3);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if 1
            % Synchrosqueezed wavelet transform (SSWT)
            
            % set parameters
            opts = struct();
            opts.motherwavelet = 'Cinfc' ;
            opts.CENTER = 1 ;
            opts.FWHM = 0.3 ;
            opts.dim = 2;%1; 5% maximum order
            
            lowfreq = 0 ;
            highfreq = 130 ;
            Hz = 2^10 ;
            L = 1 ;
            time = [1/Hz:1/Hz:L]' ;
            N = length(time) ;
            alpha = 0.2 ;
            MT = 10;
            Smooth = 1;
            Hemi = 1;
            scrsz = get(0,'ScreenSize');
            xm = x;
            
            % apply the ConceFT_CWT
            [tfrsq, ConceFT, tfrsqtic] = ConceFT_CWT(time, xm, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi) ;
            
            EMDrec2 = EMDMat(abs(ConceFT),InstFreq,highfreq);
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
            WinLen = 401;
            alpha = 0.025/Hz ;
            initstate(1) ;
            xm = x ;
            time = [1/Hz:1/Hz:L]' ;
            MT = 10; %MT = 1: ordinary SST; MT > 1: ConceFT
            hop = 1;
            dim = 2; %1; 5% maximum order
            supp = 6; %half time support (>= 6 recommended)
            Smooth = 1; % 0 or 1
            Hemi = 1; % 0 or 1
            
            % apply the ConceFT_STFT
            [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(xm,  0, 0.5*highFq/Hz*2, alpha, hop, WinLen, dim, supp, MT, Smooth, Hemi) ;
            
            EMDrec3 = EMDMat(abs(ConceFT),InstFreq,highfreq);
        end
        recVec = recVec + [EMDrec1 EMDrec2 EMDrec3]/numTest;
    end
    recMat = [recMat; recVec];
    head = sprintf('results/recQuan4_NTest_%d.mat',numTest);
    save(head,'recMat');
end
else
    head = sprintf('results/recQuan4_NTest_%d.mat',numTest);
    load(head);
end

%%
if 1
    [sm nplot] = size(recMat);
    
    pic = figure;
    hold on;
    h = zeros(1, 3);
    x = NMvec;
    x = x(:);
    y = recMat(:,1);
    h(1) = plot(x,y,'-*k');
    y = recMat(:,2);
    h(2) = plot(x,y,'-^b');
    y = recMat(:,3);
    h(3) = plot(x,y,'-om');
    hold off;
    
    legend(h,'SSWPT10','ConceFT-CWT','ConceFT-STFT','Location','southeast');
    axis square;
    hold off;
    xlabel('\sigma^2');ylabel('EMD');
    
    %     set(gca, 'FontSize', 18);
    %     b=get(gca);
    %     set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    %     str = sprintf('results/compRS%d',numTest);
    %     print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
