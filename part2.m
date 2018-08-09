clear all
%  load('price.mat')
load('dash.mat')
load('crypto.mat')
price=dash.Close;
price_7     =price(1:7:end);
price_14    =price(1:14:end);
price_30    =price(1:30:end);

% [h1,pValue1]=kpsstest(price,'alpha',0.05);
% [h7,pValue7]=kpsstest(price_7,'alpha',0.05);
% [h30,pValue30]=kpsstest(price_30,'alpha',0.05);
% 
% %% ADF test
% [h_1,pValue_1] = adftest(price,'alpha',0.05);
% [h_7,pValue_7] = adftest(price_7,'alpha',0.05);
% [h_30,pValue_30] = adftest(price_30,'alpha',0.05);


%% 

returns      = diff(price)./price(1:end-1,:);
returns7    =(price(8:end)-price(1:end-7))./price(1:end-7,:);
returns14   =(price(15:end)-price(1:end-14))./price(1:end-14,:);
returns30   =(price(31:end)-price(1:end-30))./price(1:end-30,:);
returns_7     =returns7(1:7:end);
returns_14    =returns14(1:14:end);
returns_30    =returns30(1:30:end);

lreturns      = price2ret(price);
lreturns7    =log(price(8:end))-log(price(1:end-7));
lreturns14   =log(price(15:end))-log(price(1:end-14));
lreturns30   =log(price(31:end))-log(price(1:end-30));
lreturns_7     =lreturns7(1:7:end);
lreturns_14    =lreturns14(1:14:end);
lreturns_30    =lreturns30(1:30:end);

% lreturns = lreturns_7;
% [p,value] = vratiotest(lreturns);
% [p7,value7] = vratiotest(lreturns_7);
% [p30,value30] = vratiotest(lreturns_30);

%% ARIMA 
% use daily log returns
figure
autocorr(lreturns) % MA(2)
xlabel('lag(daily)')
print('auto.eps','-depsc')
%% arima

ToEstMdl = arima('MALags',1:1,'SMALags',1);
EstMdl = estimate(ToEstMdl,lreturns);
[YF,YMSE] = forecast(EstMdl,20,'Y0',lreturns);
%% residual plot of ARIMA

Mdl = arima('Constant',0.00560328  ,'AR',-0.238106, 'variance',0.0073803 );
rng 'default';
Y = lreturns;

E = infer(Mdl,lreturns,'Y0',lreturns);
figure;
plot(E,'+');
xlabel('date')
ylabel('residuals')
title ('ARIMA (1,1,1) residual plot (Dash)');
print('residual plot.eps','-depsc')

%% forecast 
nr= lreturns(1:997)
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl1 = estimate(Mdl,lreturns);
v1 = infer(EstMdl1,lreturns);


EstMdl = estimate(Mdl,nr);

numPeriods = 427;
vF = forecast(EstMdl,numPeriods,'Y0',nr);
v = infer(EstMdl,lreturns);

figure;
plot(t(2:end),v1,'k:','LineWidth',2);
hold on;
plot(t(end)-426:t(end),vF,'r','LineWidth',2);
title('Forecasted Conditional Variances of Log-returns (Dash)');
xlim([0,length(t)]),datetick
ylabel('Conditional variances');
xlabel('Year');
legend({'Estimation sample cond. var.','Forecasted cond. var.'},...
    'Location','Best');
print('70%.eps','-depsc')

%% GARCH modeling rolling log-returns(dash)
load('t.mat')
 
nr = lreturns;

h = archtest(lreturns);
% [h,p] = lbqtest(lreturns,'Lags',1);

figure;
plot(t(1:1424),nr);
hold on;
plot([t(1) t(end-1)],[0 0],'r:'); % Plot y = 0
hold off;
title('Danish Nominal Stock Returns');
ylabel('Nominal return (%)');
xlabel('Year');

Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
% Mdl.Distribution = 't';
EstMdl = estimate(Mdl,lreturns);

numObs = numel(lreturns); % Sample size (T)
numPaths = 30;            % Number of paths to simulate
rng(1);                   % For reproducibility
[VSim,YSim] = simulate(EstMdl,numObs,'NumPaths',numPaths);
VSimBar = mean(VSim,2);
VSimCI = quantile(VSim,[0.025 0.975],2);
YSimBar = mean(YSim,2);
YSimCI = quantile(YSim,[0.025 0.975],2);

figure;
% subplot(2,1,1);
h1 = plot(t(2:1425),VSim,'Color',0.8*ones(1,3));
hold on;
h2 = plot(t(2:1425),VSimBar,'k--','LineWidth',2);
h3 = plot(t(2:1425),VSimCI,'r--','LineWidth',2);
hold off;
xlim([0,length(t)]),datetick
title('Simulated Conditional Variances');
ylabel('Cond. var.');
xlabel('Year');
legend([h1(1) h2 h3(1)],{'Simulated path' 'Variance' 'Confidence bounds'},...
    'FontSize',7,'Location','NorthWest');
print('gggg.eps','-depsc')


%% 

% figure;
% % subplot(2,1,2);
% h1 = plot(t(2:1425),YSim,'Color',0.8*ones(1,3));
% hold on;
% h2 = plot(t(2:1425),YSimBar,'k--','LineWidth',2);
% h3 = plot(t(2:1425),YSimCI,'r--','LineWidth',2);
% hold off;
% xlim([0,length(t)]),datetick
% title('Simulated Nominal Returns');
% ylabel('Nominal return (%)');
% xlabel('Year');
% legend([h1(1) h2 h3(1)],{'Simulated path' 'Mean' 'Confidence bounds'},...
%     'FontSize',7,'Location','NorthWest');

%% 
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl1 = estimate(Mdl,lreturns);
v1 = infer(EstMdl1,lreturns);
vF=zeros(427,1);
for i=997:1423
   nr= lreturns(1:i);
   EstMdl2 = estimate(Mdl,nr,'Display','off');
   vF(i-996,1) = forecast(EstMdl2,1,'Y0',nr);
end

%% 
figure;
hold on 
%  plot(t(1:1424),lreturns,'k');
plot(t(1:1424),lreturns,'Color',0.8*ones(1,3));
plot(t(2:end),v1,'k:','LineWidth',2);
xlabel('Year');
legend({'log-returns','Estimation sample cond. var.'},...
    'Location','Best');
title('Compare Conditional Variances with log-returns');
print('Compare Conditional Variances with log-returns.eps','-depsc')

%% 

figure;
hold on;
plot(t(2:end),v1,'k:','LineWidth',2);
plot(t(end)-426:t(end),vF,'r','LineWidth',1);
plot([t(1) t(end-1)],[0 0],'r:'); % Plot y = 0
title('Forecasted Conditional Variances of log-returns');
xlim([0,length(t)]),datetick
ylabel('Conditional variances');
xlabel('Year');
legend({'Estimation sample cond. var.','Forecasted cond. var.'},...
    'Location','Best');
print('Forecasted Conditional Variances of log-returns.eps','-depsc')

%% Hurst plot

t=1:4196;
K_UP=zeros(15,1);
K_DOWN=zeros(15,1);
K1=zeros(15,1);
    for tao=1:15
        for t=1:4196-tao
        K_up=abs((price(t+tao,:))-(price(t,:)));
        K_UP(tao,:)=K_UP(tao,:)+K_up;
        K_down=abs(price(t,:));
        K_DOWN(tao,:)=K_DOWN(tao,:)+K_down;
        end
        K1(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
    end
  log_tao=log(1:15);
  log_K1=log(K1);
  
 %q=2
 K_UP=zeros(15,1);
 K_DOWN=zeros(15,1);
 K2=zeros(15,1);
    for tao=1:15
        for t=1:4196-tao
        K_up(tao)=abs(price(t+tao,:)-price(t,:)).^2;
        K_UP(tao,:)=K_UP(tao,:)+K_up(tao);
        K_down(tao)=abs(price(t,:)).^2;
        K_DOWN(tao,:)=K_DOWN(tao,:)+K_down(tao);
        end
        K2(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
    end
  log_tao=log(1:15);
  log_K2=log(K2);
  
 %q=3
K_UP=zeros(15,1);
K_DOWN=zeros(15,1);
K3 =zeros(15,1);
    for tao=1:15
        for t=1:4196-tao
        K_up(tao)=abs(price(t+tao,:)-price(t,:)).^3;
        K_UP(tao,:)=K_UP(tao,:)+K_up(tao);
        K_down(tao)=abs(price(t,:)).^3;
        K_DOWN(tao,:)=K_DOWN(tao,:)+K_down(tao);
        end
        K3(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
    end
  log_tao=log(1:15);
  log_K3=log(K3);
  
 %plot GHE
figure
% plot(log_tao,log_K1,'-om',log_tao,log_K2,'-or',log_tao,log_K3,'-ob')
plot(log_tao,log_K1,'-o',log_tao,log_K2,'-o',log_tao,log_K3,'-o')

% xlim([0,3]);
  ylim([-15,0]);
legend('q=1','q=2','q=3')
title('Hurst exponent(D3446)')
xlabel('log(tao)')
ylabel('log(kq(tao))')
print('hurst(d3446).eps','-depsc')
%% Hurst exponent

hh1=GenHurst(price,1);
hh2=GenHurst(price,2);
hh3=GenHurst(price,3);

hh =zeros(30,1);

for p=0:1:30
    hh(p+1,:)=p*GenHurst(price,p/30);
end

p=0:0.1:3;
figure
plot(p,hh,'-o')
title('Hurst exponent(Dash)')
xlabel('q')
ylabel('q*H(q)')
print('hurst(Dash).eps','-depsc')

%% KPSS test
[h1,pValue1]=kpsstest(lreturns,'alpha',0.05);
[h7,pValue7]=kpsstest(lreturns_7,'alpha',0.05);
[h30,pValue30]=kpsstest(lreturns_30,'alpha',0.05);

%% ADF test
[h_1,pValue_1] = adftest(lreturns,'alpha',0.05);
[h_7,pValue_7] = adftest(lreturns_7,'alpha',0.05);
[h_30,pValue_30] = adftest(lreturns_30,'alpha',0.05);
%% pacf
figure
parcorr(lreturns)

%% acf daily
figure
subplot(3,1,1)
autocorr(lreturns)
title('KO autocorrelation of log-returns','fontsize',14)
xlabel('lags (days)')
% acf=autocorr(lreturns);
subplot(3,1,2)
autocorr(abs(lreturns))
title('KO autocorrelation of |log-returns|','fontsize',14)
xlabel('lags (days)')
subplot(3,1,3)
autocorr(lreturns.^2)
title('KO autocorrelation of log-returns^2','fontsize',14)
 xlabel('lags (days)')
% ylabel('autocorrelation','fontsize',14)
% title('KO autocorrelation of log-returns(daily)','fontsize',14)
% set(gca,'fontsize',14)
% legend({'log-returns','|log-returns|','log-returns^2'})
print('acf.eps','-depsc')

%% acfweekly
figure
subplot(3,1,1)
autocorr(lreturns_7)
title('KO autocorrelation of log-returns','fontsize',14)
xlabel('lags (weeks)')
% acf=autocorr(lreturns);
subplot(3,1,2)
autocorr(abs(lreturns_7))
title('KO autocorrelation of |log-returns|','fontsize',14)
xlabel('lags (weeks)')
subplot(3,1,3)
autocorr(lreturns_7.^2)
title('KO autocorrelation of log-returns^2','fontsize',14)
 xlabel('lags (weeks)')
% ylabel('autocorrelation','fontsize',14)
% title('KO autocorrelation of log-returns(daily)','fontsize',14)
% set(gca,'fontsize',14)
% legend({'log-returns','|log-returns|','log-returns^2'})
print('acfweek.eps','-depsc')

%% acfmonthly
figure
subplot(3,1,1)
autocorr(lreturns_30)
title('KO autocorrelation of log-returns','fontsize',14)
xlabel('lags (months)')
% acf=autocorr(lreturns);
subplot(3,1,2)
autocorr(abs(lreturns_30))
title('KO autocorrelation of |log-returns|','fontsize',14)
xlabel('lags (months)')
subplot(3,1,3)
autocorr(lreturns_30.^2)
title('KO autocorrelation of log-returns^2','fontsize',14)
 xlabel('lags (months)')
% ylabel('autocorrelation','fontsize',14)
% title('KO autocorrelation of log-returns(daily)','fontsize',14)
% set(gca,'fontsize',14)
% legend({'log-returns','|log-returns|','log-returns^2'})
print('acfmonth.eps','-depsc')

%% 
figure
subplot(2,1,1)
[a1,lags1] = autocorr(lreturns(:,1),250);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns(:,1)),250);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns(:,1).^2,250);
plot(lags3,a3,'-b')
plot([0 300],[0 0])
axis([0 250 -.15 1])
xlabel('lags (days)')
ylabel('autocorrelation')
title('KO autocorrelation of log-returns(daily)')
legend({'log-returns','|log-returns|','log-returns^2'})

subplot(2,1,2)
[a1,lags1] = autocorr(lreturns_7(:,1),200);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns_7(:,1)),200);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns_7(:,1).^2,200);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 200 -.15 1])
xlabel('lags (weeks)')
ylabel('autocorrelation')
title('KO autocorrelation of log-returns(weekly)')
legend({'log-returns','|log-returns|','log-returns^2'})
print('autocorrelationDay.eps','-depsc')

figure
subplot(2,1,1)
[a1,lags1] = autocorr(lreturns_30(:,1),40);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns_30(:,1)),40);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns_30(:,1).^2,40);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 40 -.15 1])
xlabel('lags (months)')
ylabel('autocorrelation')
title('KO autocorrelation of log-returns(monthly)')
legend({'log-returns','|log-returns|','log-returns^2'})
print('autocorrelation30Day.eps','-depsc')

%% autocorrelogram
figure
[a1,lags1] = autocorr(lreturns(:,1),250);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns(:,1)),250);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns(:,1).^2,250);
plot(lags3,a3,'-b')
plot([0 300],[0 0])
axis([0 250 -.15 1])
xlabel('lags (days)','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('KO autocorrelation of log-returns(daily)','fontsize',14)
set(gca,'fontsize',14)
legend({'log-returns','|log-returns|','log-returns^2'})
print('autocorrelation1Day.eps','-depsc')

%% 

figure
[a1,lags1] = autocorr(lreturns_7(:,1),200);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns_7(:,1)),200);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns_7(:,1).^2,200);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 200 -.15 1])
xlabel('lags (weeks)','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('KO autocorrelation of log-returns(weekly)','fontsize',14)
set(gca,'fontsize',14)
legend({'log-returns','|log-returns|','log-returns^2'})
print('autocorrelation7Day.eps','-depsc')

figure
[a1,lags1] = autocorr(lreturns_30(:,1),40);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(lreturns_30(:,1)),40);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(lreturns_30(:,1).^2,40);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 40 -.15 1])
xlabel('lags (months)','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('KO autocorrelation of log-returns(monthly)','fontsize',14)
set(gca,'fontsize',14)
legend({'log-returns','|log-returns|','log-returns^2'})
print('autocorrelation30Day.eps','-depsc')

%% scale & rescale

t1=find(diff(price)~=0,1,'first')+1;
t2=find(diff(price)~=0,1,'last');

pr=price(t1:t2); 

mm=[1 7 14 30];
for m=mm,lr{m} = [];end

for m=mm
    lr{m} = [lr{m};log(pr((m+1):m:end))-log(pr(1:m:(end-m)))];
end

h=GenHurst(pr,2);
bins=[];
b1=[];
bins1=[];
freqs=[];
for m=mm
    s(m) = std(lr{m});
    [f,b]=hist(lr{m},30);
    b1=[b1,b'];
    bins = [bins,b'/m];
    bins1 = [bins1,b'/m^h];
    freqs = [freqs,f'/max(f)];
end
%% 

figure
semilogy(bins,freqs,'-+')
%text(0.002,0.2,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.5 +0.5 0 1.1])
ylim([10^-4,10^0]);
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('Scale(D3446)')
print('Scale(D3446).eps','-depsc')
%% 

figure
semilogy(b1,freqs,'-+')
%text(0.002,0.2,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.23 +0.23 0 1.1])
ylim([10^-4,10^0]);
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('Scale(D3446)')
print('Scale(D3446).eps','-depsc')
%% 

figure
semilogy(bins1,freqs,'-+')
hold on
x=-1:0.001:1;
% nor=fitdist(lr{1},'normal');
normdist = fitdist(lreturns,'normal');
stamdist = fitdist(lreturns,'tlocationscale');
y1=pdf(normdist,x);
y2=pdf(stamdist,x);
semilogy(x,y2/max(y2)*0.55,'--k')
text(-0.1,10^-1,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.13 +0.13 5e-4 1.1])
ylim([10^-4,10^0]);
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return/¦Ó^H(2)','fontsize',14)
ylabel('relative freq','fontsize',14)
title('RetDistrScaled(D3446)')
print('RetDistrScaled(D3446).eps','-depsc')


%% scale 
% 
% t1=find(diff(price)~=0,1,'first')+1;
% t2=find(diff(price)~=0,1,'last');
% 
% pr=price(t1:t2); 
% 
% syb = 'v<>^s';
% col = 'rbkgm';
% k=0;
% mm=[1 7 14 30];
% 
% figure
% for m=mm
%     k=k+1;
%     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
%     [f,b]=hist(lr,30);
%     semilogy(b,f/max(f),['-',syb(k),col(k)])
%     hold on
% end
% 
% axis([-0.5 +1 1e-4 1.1])
% legend(num2str(mm'))
% set(gca,'fontsize',14)
% xlabel('lor-return','fontsize',14)
% ylabel('relative freq','fontsize',14)
% title('AAPL 1 sec data 3 Dic 2012')
% print('1secRetDistr.eps','-depsc')
% 
% %% rescale
% h=GenHurst(pr,2);
% k=0;
% figure
% for m=mm
%     k=k+1;
%     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
%     [f,b]=hist(lr,30);
%     semilogy(b/m.^h,f/max(f),['-',syb(k),col(k)])
%     hold on
% end
% % text(0.0005,0.05,['H(2) = ',num2str(h)],'fontsize',14)
% axis([-0.5 +0.5 1e-4 1.1])
% legend(num2str(mm'))
% set(gca,'fontsize',14)
% xlabel('log-return','fontsize',14)
% ylabel('relative freq','fontsize',14)
% title('AAPL 1 sec data 3 Dic 2012')
% print('1secRetDistrScaled.eps','-depsc')

%% complementary cumulative distribution
figure,
lreturns1 = lreturns(:,1);
lreturns7 = lreturns_7(:,1);

loglog(sort(lreturns1),1-[1:(length(lreturns1))]...
    /length(lreturns1),'+b')
hold on
loglog(sort(lreturns7),1-[1:(length(lreturns7))]...
    /length(lreturns7),'xr')

x = [max(lreturns1)/1000:max(lreturns1)/1000:max(lreturns1)];
axis([1e-5 0.5 1e-4 1])
legend({'pos ret','neg ret'})
title('complemetary cumulative log-return distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
print('fat_tailsCum.eps','-depsc')

%% alpha and hurst relationship

H1=zeros(25,1);
qH1=zeros(25,1);
i=1;
alpha=3.6;
for q=0:0.25:6 
    
    H1(i,1)=GenHurst(cumsum(lreturns),q);
   
    qH1(i,1)=H1(i,1)*q;
    i=i+1;
end

figure,
subplot(2,1,1)
err=0.05*qH1;
errorbar([0:0.25:6],qH1,err)
hold on

x1=3.6:0.01:6;
y1=ones(length(x1),1)*3.6/2;
plot(x1,y1,'linewidth',1.5)

y2=0:0.01:2.5;
x2=ones(length(y2),1)*alpha
plot(x2,y2,'linewidth',1.5)

x3=0:0.01:3.6;
y3=0.5*x3;
plot(x3,y3,'linewidth',1.5)

xlabel('q')
ylabel('qH(q)')
text(1,0.3,'<|r(\tau)|^q>~\tau^{-qH}','fontsize',18)
text(4.3,1.3,'<|r(\tau)|^q>~\tau','fontsize',18)
text(2.3,2, '\alpha=3.7174','fontsize',18)
title('Relationship between hurst exponent and tail exponent (D3446)','fontsize',14)
print('d3446 Relationship between hurst exponent and tail exponent.eps','-depsc')