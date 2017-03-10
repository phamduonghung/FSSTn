function fig2()

% fig2 : draws Figure 2 of paper "High-Order Synchrosqueezing Transform for
% Multicomponent Signals Analysis - With an Application to Gravitational-Wave Signal, by PHAM and Meignen.

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Journal_IEEE_2016/Tex/Figures';

% Parameters
gamma = 0;
N = 1024;

t = (0:N-1)/N;f = 0:N/2-1;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;
NMvec = -5:5:5;
sigmaR = 0.01:0.01:0.2;
Ren = zeros(length(sigmaR),length(NMvec));
sigma_optR = zeros(1,length(NMvec));

%% Test signal 
[a1,a2,if1,if2,s1,s2,s] = signal_test(t);

%% Iteration for different sigma valeurs

% No noise
Ren_no = zeros(length(sigmaR),1);
for j = 1:length(sigmaR)
% SST
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s,gamma,sigmaR(j),ft,bt);
Ren_no(j) =renyi_entropy(abs(STFT),t,f',3);
end
sigma_optR_no = sigmaR(Ren_no==min(Ren_no))

%Noise
for cntNM = 1:length(NMvec) 
    % loop on noise
    
    sg = sqrt( var(s)*10^(-NMvec(cntNM)/10))
    % set noise
     b = sg*randn(size(s));
     sb = s+b;
     SNR = 10*log10((var(s(:)))/var(b(:)))
     
    for j = 1:length(sigmaR)
        % SST
        [STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(sb,gamma,sigmaR(j),ft,bt);
        Ren(j,cntNM) =renyi_entropy(abs(STFT),t,f',3);
    end
        sigma_optR(cntNM) = sigmaR(Ren(:,cntNM)==min(Ren(:,cntNM)))
end


%% Display the concentration measure!!

FigHandle113 = figure; 
set(FigHandle113, 'Position', [350,350, 1100, 900]);

    plot(sigmaR,Ren_no,'m-x',sigmaR,Ren(:,1),'k*-',sigmaR,Ren(:,2),'b-s',sigmaR,Ren(:,3),'r-o'); 
    
    strmin = [' \leftarrow H_{R,min}=' num2str(min(Ren_no),'%.2f') ' & \sigma_{opt} = ' num2str(min(sigma_optR_no))];
    text(sigma_optR_no,min(Ren_no)-0.1,strmin,'HorizontalAlignment','left','color','m');
    
    
    strmin = [' \leftarrow H_{R,min}=' num2str(min(Ren(:,1)),'%.2f') ' & \sigma_{opt} = ' num2str(min(sigma_optR(1)))];
    text(sigma_optR(1),min(Ren(:,1))-0.1,strmin,'HorizontalAlignment','left','color','k');
    
    strmin = [' \leftarrow H_{R,min}=' num2str(min(Ren(:,2)),'%.2f') ' & \sigma_{opt} = ' num2str(min(sigma_optR(2)))];
    text(sigma_optR(2),min(Ren(:,2))-0.1,strmin,'HorizontalAlignment','left','color','b');
    
    
    strmin = [' \leftarrow H_{R,min}=' num2str(min(Ren(:,3)),'%.2f') ' & \sigma_{opt} = ' num2str(min(sigma_optR(3)))];
    text(sigma_optR(3),min(Ren(:,3))-0.1,strmin,'HorizontalAlignment','left','color','r');
    
    ylim([min(Ren_no)-0.4, max(Ren(:,3))+1.2]);
    xlim([min(sigmaR), max(sigmaR)]);
    legend('noise-free','-5 dB','0 dB','5bB','Location','northeast');

 ylabel('Renyi entropy');
 xlabel('\sigma');
 set(findall(0,'type','Line'),'LineWidth',1);
 set(findall(FigHandle113,'-property','FontSize'),'FontSize',24)

%%%%%%%%%%%%%%%%%%%%%% print signal
export_fig(FigHandle113, ... % figure handle
    sprintf('%s/Renyi_entropy', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
 

end
