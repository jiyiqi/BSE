clear
close all
%%%%%%%%%%%%%%%
% load daily prices for Apple 
load('AAPL_daily_1990_2012.mat') % data from http://www.dukascopy.com/swiss/english/marketwatch/historical/
%%%%%%%%%%%%%%%%
% extract trading days only
t = find(Volume>1e11);
Day = Day(t);
Month = Month(t);
Year = Year(t);
Dates = Dates(t);
ClosingPrice = ClosingPrice(t);
Volume = Volume(t);
%%%%%%%%%%%%%%%
% plot prices
figure
plot(Dates,ClosingPrice),datetick
title('AAPL prices')
xlabel('date','fontsize',14)
ylabel('price','fontsize',14)
set(gca,'fontsize',14)
print('prices.eps','-depsc')

%%%%%%%%%%%%%%%
% plot volumes
figure
plot(Dates,Volume),datetick
title('AAPL volumes')
xlabel('date','fontsize',14)
ylabel('volume','fontsize',14)
set(gca,'fontsize',14)
print('volumes.eps','-depsc')

%%%%%%%%%%%%%%%
% compute daily returns
returns = diff(ClosingPrice);
%%%%%%%%%%%%%%%
% compute daily log-returns
lreturns = diff(log(ClosingPrice));
%%%%%%%%%%%%%%%
% plot returns
figure
plot(Dates(2:end),lreturns),datetick
title('AAPL log-returns')
xlabel('date','fontsize',14)
ylabel('log-return','fontsize',14)
set(gca,'fontsize',14)
print('log_returns.eps','-depsc')


%%%%%%%%%%%%%%%%%%%
% stylized fact -1-
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% compute and plot return distribution
[freq,bin]=hist(lreturns,200);
figure
bar(bin,freq/sum(freq))
%plot(bin,freq/sum(freq),'o')
hold on
% compare with normal distribution
m = mean(lreturns);
s = std(lreturns);
x = [min(lreturns):(max(lreturns)-min(lreturns))/1000:max(lreturns)];
plot(x,normpdf(x,m,s)*(bin(2)-bin(1)),'-m','linewidth',2)
axis([-0.2 0.2 0 max(freq/sum(freq))*1.1])
legend({'log-return','normal'})
title('AAPL log-return relative frequency distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)
set(gca,'fontsize',14)
print('fat_tails.eps','-depsc')
%%%%%%%%%%%%%%%
% compute and plot complemetary cumulative return distribution (rank-frequency plot method)
figure,
loglog(sort(lreturns(lreturns>0)),1-[1:(length(lreturns(lreturns>0)))]/length(lreturns(lreturns>0)),'+b')
hold on
loglog(sort(-lreturns(lreturns<0)),1-[1:(length(lreturns(lreturns<0)))]/length(lreturns(lreturns<0)),'xr')
% compare with normal distribution
x = [max(lreturns)/1000:max(lreturns)/1000:max(lreturns)];
loglog(x,1-(normcdf(x,m,s)-0.5)*2,'-m','linewidth',2)
axis([1e-5 0.5 1e-4 1])
legend({'pos ret','neg ret','normal'})
title('AAPL complemetary cumulative log-return distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
print('fat_tailsCum.eps','-depsc')



%%%%%%%%%%%%%%%%%%%
% stylized fact -2-
%%%%%%%%%%%%%%%%%%%
figure
[a1,lags1] = autocorr(returns,250);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(returns),250);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(returns.^2,250);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 250 -.15 1])
xlabel('lags (days)','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('AAPL autocorrelation of returns','fontsize',14)
set(gca,'fontsize',14)
legend({'returns','|returns|','returns^2'})
print('autocorrelation1Day.eps','-depsc')

%%%%%%%%%%%%%%%%%%%
% stylized fact -3-
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% compute volatility over 20 days using the cumsum trick
m = 20;
%x = (returns./price(2:end));
x  = returns;
Y1 = cumsum(x,1) ;
Y2 = cumsum(x.^2,1) ;
volatility = sqrt((Y2((m+1):end)-Y2(1:(end-m)))/m-((Y1((m+1):end)-Y1(1:(end-m)))/m).^2) ;
x  = Volume(2:end);
Y1 = cumsum(x,1) ;
Y2 = cumsum(x.^2,1) ;
meanVolume = sqrt((Y2((m+1):end)-Y2(1:(end-m)))/m-((Y1((m+1):end)-Y1(1:(end-m)))/m).^2) ;
figure
subplot(2,1,1)
plot(Dates((m+2):end),volatility),datetick
xlabel('date','fontsize',14)
ylabel('volatility','fontsize',14)
set(gca,'fontsize',14)
subplot(2,1,2)
plot(Dates((m+2):end),meanVolume),datetick
xlabel('date','fontsize',14)
ylabel('volume','fontsize',14)
set(gca,'fontsize',14)
print('volatility.eps','-depsc')
%%
c=corrcoef(meanVolume,volatility)
figure
plot((meanVolume),(volatility),'+b')
text(0.6*max(meanVolume),0.6*max(volatility),['correlation coeff =',num2str(c(1,2))],'fontsize',14)
xlabel('volume','fontsize',14)
ylabel('volatility','fontsize',14)
set(gca,'fontsize',14)
print('volume_volatility.eps','-depsc')
%%

Dates0= Dates;
ClosingPrice0=ClosingPrice;
load MSFT_Daily_1990-2012.mat
match = zeros(1,length(Dates));
for t=1:length(Dates0)
    match = match+(floor(Dates) == floor(Dates0(t)));
end

c=corrcoef(diff(log(ClosingPrice0)),diff(log(ClosingPrice(match>0))));
figure
plot(Dates0,ClosingPrice0,'-r')
hold on
plot(Dates(match>0),ClosingPrice(match>0),'-b')
text(min(Dates0)*1.001,0.6*max(ClosingPrice0),['correlation coeff =',num2str(c(1,2))],'fontsize',14)
datetick
title('AAPL - MSFT','fontsize',14)
xlabel('date','fontsize',14)
ylabel('prices','fontsize',14)
set(gca,'fontsize',14)
print('correlations.eps','-depsc')

