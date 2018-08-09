%% load the data
 clear all, close all
 warning off all
 rng('default');
 load('dash.mat')
 load('t.mat')
 Volume1 = dash.Volume;
 price_series1=dash.Close;

% Compute Daily Returns
returns      = diff(price_series1)./price_series1(1:end-1,:);
%´ÓÎ²µ½Í·
size_returns = length(returns);

% Compute Daily log-returns
lreturns      = price2ret(price_series1);
size_lreturns = length(lreturns);


week = 1:7:length(t);
price_seriesw = price_series1(week);
returns_7   = (price_seriesw(2:end)-price_seriesw(1:end-1))./price_seriesw(1:end-1);
size_returns_7 = length(returns_7);
lreturns_7  = log(price_seriesw(2:end))-log(price_seriesw(1:end-1));
%% sample statistics
display('------------------------------')
display('Sample Statistics:')
mean_lr     = mean(lreturns(:,1));
display(sprintf('mean_lr  : %.9f', mean_lr))
std_lr      = std(lreturns(:,1));
display(sprintf('std_lr   : %.9f', std_lr))
min_lr      = min(lreturns(:,1));
display(sprintf('min_lr   : %.9f', min_lr))
max_lr      = max(lreturns(:,1));
display(sprintf('max_lr   : %.9f', max_lr))
med_lr       = median(lreturns(:,1));
display(sprintf('median_lr   : %i', med_lr))
ske_lr      = skewness(lreturns(:,1));
display(sprintf('ske_lr   : %.9f', ske_lr))
kurt_lr     = kurtosis(lreturns(:,1));
display(sprintf('kurt_lr  : %.9f', kurt_lr))
ad_lr       = adtest(lreturns(:,1));
display(sprintf('ad_lr   : %i', ad_lr))
display('------------------------------')
kst_lr      = kstest(lreturns(:,1));
display(sprintf('kst_lr   : %i', kst_lr))
display('------------------------------')

m = mean(lreturns(:,1));
s = std(lreturns(:,1));
x = [min(lreturns(:,1)):(max(lreturns(:,1))-min(lreturns(:,1)))...
     /1000:max(lreturns(:,1))];
 
pd=fitdist(lreturns,'stable');
%% comparison with T Location-Scale Distribution
figure;

h1=histfit(lreturns,500,'stable');
set(h1(2),'Color','r');
hold on
h2=histfit(lreturns,500,'tlocationscale');
set(h2(1),'Visible','off');
set(h2(2),'Color','g');
legend({'log-return frequencey','Student t'})
title('log-return relative frequency distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)
set(gca,'fontsize',14)
hold on
figure
histfit(lreturns,500,'normal')
legend('Log-Returns frequency','Gaussian')
title('Normal fitted to daily log returns')
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)

% %% comparison with Gaussian
% 
% [freq,bin]=hist(lreturns(:,1),500);
% figure
% bar(bin,freq/sum(freq)) 
% hold on
% plot(x,normpdf(x,m,s)*(bin(2)-bin(1)),'-m','linewidth',2)
% axis([-0.2 0.2 0 max(freq/sum(freq))*1.1])
% legend({'log-return','normal'})
% title('log-return relative frequency distribution','fontsize',14)
% xlabel('log-return','fontsize',14)
% ylabel('relative frequency','fontsize',14)
% set(gca,'fontsize',14)

%% complementary stable cumulative distribution
figure,
lreturns1 = lreturns(:,1);
loglog(sort(lreturns1(lreturns1>0)),1-[1:(length(lreturns1(lreturns1>0)))]...
    /length(lreturns1(lreturns1>0)),'+b')
hold on
loglog(sort(-lreturns1(lreturns1<0)),1-[1:(length(lreturns1(lreturns1<0)))]...
    /length(lreturns1(lreturns1<0)),'xr')

x = [max(lreturns1)/1000:max(lreturns1)/1000:max(lreturns1)];
loglog(x,1-(cdf(pd,x)-0.5)*2,'-m','linewidth',2)
axis([1e-5 0.5 1e-4 1])
legend({'pos ret','neg ret','stable distribution'},'Location','northwest')
title('complemetary cumulative log-return distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
print('ccdf (stable).eps','-depsc')

%% complementary Normal cumulative distribution
figure,
lreturns1 = lreturns(:,1);
loglog(sort(lreturns1(lreturns1>0)),1-[1:(length(lreturns1(lreturns1>0)))]...
    /length(lreturns1(lreturns1>0)),'+b')
hold on
loglog(sort(-lreturns1(lreturns1<0)),1-[1:(length(lreturns1(lreturns1<0)))]...
    /length(lreturns1(lreturns1<0)),'xr')
x = [max(lreturns1)/1000:max(lreturns1)/1000:max(lreturns1)];
loglog(x,1-(normcdf(x,m,s)-0.5)*2,'-m','linewidth',2)
axis([1e-5 0.5 1e-4 1])
legend({'pos ret','neg ret','normal'})
title('complemetary cumulative log-return distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
print('fat_tailsCumnormal.eps','-depsc')

%% qq plot
figure,
qqplot(lreturns(:,1)),

figure,
PD = fitdist(lreturns(:,1),'GeneralizedExtremeValue'); 
qqplot(lreturns(:,1),PD)

figure,
PD = fitdist(lreturns(:,1),'tLocationScale');
qqplot(lreturns(:,1),PD)

