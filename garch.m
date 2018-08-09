clear all
load('price.mat')
 
price_7     =price(1:7:end);
price_14    =price(1:14:end);
price_30    =price(1:30:end);

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
t=1:4196;
%% ARIMA 
% use daily log returns
figure
autocorr(lreturns) % MA(2)
xlabel('lag(daily)')
print('auto3446.eps','-depsc')
%% arima

ToEstMdl = arima('MALags',2:1,'SMALags',2);
EstMdl = estimate(ToEstMdl,lreturns);
[YF,YMSE] = forecast(EstMdl,20,'Y0',lreturns);
%% residual plot of ARIMA

Mdl = arima('Constant',0.000269271  ,'AR',-0.238106, 'variance',0.000141002 );
rng 'default';
Y = lreturns(400:2000);

E = infer(Mdl,lreturns(400:2000),'Y0',lreturns(400:2000));
figure;
plot(E,'+');
xlabel('date')
ylabel('residuals')
title ('ARIMA (2,1,2) residual plot (D3446)');
print('residual plot3446.eps','-depsc')

%% 

nr= lreturns(1:2938)
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl1 = estimate(Mdl,lreturns);
v1 = infer(EstMdl1,lreturns);


EstMdl = estimate(Mdl,nr);

numPeriods = 1258;
vF = forecast(EstMdl,numPeriods,'Y0',nr);
v = infer(EstMdl,lreturns);

figure;
plot(t(1:end),v1,'k:','LineWidth',2);
hold on;
plot(t(end)-1257:t(end),vF,'r','LineWidth',2);
title('Forecasted Conditional Variances of Log-returns (D3446)');
xlim([0,length(t)]),datetick
ylabel('Conditional variances');
xlabel('Year');
legend({'Estimation sample cond. var.','Forecasted cond. var.'},...
    'Location','Best');
print('70%.eps','-depsc')

%% GARCH modeling rolling log-returns(dash)
 
nr = lreturns;

h = archtest(lreturns);
% [h,p] = lbqtest(lreturns,'Lags',1);

figure;
plot(t(1:4196),nr);
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
h1 = plot(t(1:4196),VSim,'Color',0.8*ones(1,3));
hold on;
h2 = plot(t(1:4196),VSimBar,'k--','LineWidth',2);
h3 = plot(t(1:4196),VSimCI,'r--','LineWidth',2);
hold off;
xlim([0,length(t)]),datetick

title('Simulated Conditional Variances');
ylabel('Cond. var.');
xlabel('Year');
legend([h1(1) h2 h3(1)],{'Simulated path' 'Variance' 'Confidence bounds'},...
    'FontSize',7,'Location','NorthWest');
print('gg.eps','-depsc')

%% 
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl1 = estimate(Mdl,lreturns);
v1 = infer(EstMdl1,lreturns);
vF=zeros(1258,1);
for i=2938:4196
   nr= lreturns(1:i);
   EstMdl2 = estimate(Mdl,nr,'Display','off');
   vF(i-2937,1) = forecast(EstMdl2,1,'Y0',nr);
end

%% 
figure;
hold on 
%  plot(t(1:1424),lreturns,'k');
plot(t(1:4196),lreturns,'Color',0.8*ones(1,3));
plot(t(1:end),v1*30,'k:','LineWidth',2);
xlabel('Year');
legend({'log-returns','Estimation sample cond. var.'},...
    'Location','Best');
title('Compare Conditional Variances with log-returns');
print('Compare Conditional Variances with log-returns.eps','-depsc')

%% 

figure;
hold on;
plot(t(1:end),v1,'k:','LineWidth',2);
plot(t(end)-1258:t(end),vF,'r','LineWidth',1);
plot([t(1) t(end-1)],[0 0],'r:'); % Plot y = 0
title('Forecasted Conditional Variances of log-returns');
xlim([0,length(t)]),datetick
ylabel('Conditional variances');
xlabel('Year');
legend({'Estimation sample cond. var.','Forecasted cond. var.'},...
    'Location','Best');
print('Forecasted Conditional Variances of log-returns.eps','-depsc')
