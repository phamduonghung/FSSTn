function plotsparse_new(R,TvO,Tv,T,nc,es)
% plotsparse_new: plots energy coefficients (see fig2_VSST)

N = size(T,2);
fsize=15; % font size
    
% Sort
R = sort(abs(R(:)),'descend');
%T1 = T;
T = sort(abs(T(:)),'descend');
%T2=T;
%figure()
%imagesc(ones(size(T1)).*(abs(T1)>=T2(N)))
%pause
%sum(abs(T2(1:N)).^2)
%sum(sum(abs(T1).^2))
%sum(T.^2)
%plot(abs(T1(:,150)))
%pause;
Tv = sort(abs(Tv(:)),'descend');
TvO = sort(abs(TvO(:)),'descend');

% cumulative sum
R = cumsum(R.^2)/sum(R.^2);
T = cumsum(T.^2)/sum(T.^2);
%T(N)
Tv = cumsum(Tv.^2)/sum(Tv.^2);
TvO = cumsum(TvO.^2)/sum(TvO.^2);

% abscissae;
P = 0.5:0.5:nc;
P = round(P*N);
R = R(P);
T = T(P);
Tv = Tv(P);
TvO = TvO(P);
x = P/N;%+10/N;

figure();
hold on;

i = 1:length(P);
plot(x(i),R(i),'bx-','markersize',6);
plot(x(i),TvO(i),'go-','markersize',6);
plot(x(i),Tv(i),'m*-','markersize',6);
plot(x(i),T(i),'rd-','markersize',6);

legend('RM','FSST2','FSST3','FSST4', 'Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
set(gca,'fontsize',fsize,'linewidth',1)
box on;
end