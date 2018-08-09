clear all
 load('price.mat')
%  load('dash.mat')
%  price=dash.Close;


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

[alpha1, xmin1, L1]=plfit(lreturns);
[alpha2, xmin2, L1]=plfit(lreturns_7);
[alpha3, xmin3, L1]=plfit(lreturns_30);

%% Sample Statistics
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

%% comparison with Gaussian 方法1
[freq,bin]=hist(lreturns(:,1),500); 
figure
bar(bin,freq/sum(freq)) 
hold on 
plot(x,normpdf(x,m,s)*(bin(2)-bin(1)),'-m','linewidth',2) 
axis([-1 1 0 max(freq/sum(freq))*1.1]) 
xlim([-1 1.5]);
legend({'weekly log-return','normal'}) 
title('log-return relative frequency distribution','fontsize',14) 
xlabel('log-return','fontsize',14) 
ylabel('relative frequency','fontsize',14) 
set(gca,'fontsize',14)
print('comparison1 with tlocation.eps','-depsc')

figure 
bar(bin,freq/sum(freq)) 
hold on;
sdist = fitdist(lreturns,'tlocationscale');
normdist = fitdist(lreturns,'normal');
plot(x,pdf(sdist,x)*(bin(2)-bin(1)),'-g','linewidth',2);
plot(x,pdf(normdist,x)*(bin(2)-bin(1)),'-m','linewidth',2);
axis([-1 1 0 max(freq/sum(freq))*1.1]) 
xlim([-1 1.5]);
xlabel('log-return','fontsize',14) 
ylabel('relative frequency','fontsize',14) 
legend({'wwekly log-return','tlocation scale distribution','normal'})
print('comparison with tlocation.eps','-depsc')

%% comparison with normal Distribution方法2
figure
histfit(lreturns,500,'normal')
legend('Log-Returns frequency','Gaussian')
title('Normal fitted to daily log returns')
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)

 %% comparison with T Location-Scale Distribution

figure;
h1=histfit(lreturns,500,'tlocationscale');
set(h1(2),'Color','r');
hold on
h2=histfit(lreturns,500,'normal');
set(h2(1),'Visible','off');
set(h2(2),'Color','g');
legend({'log-return frequencey','T locationscale'})
title('log-return relative frequency distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)
set(gca,'fontsize',14)

%% comparison stable with t

figure;
h1=histfit(lreturns,500,'tlocationscale');
set(h1(2),'Color','r');
hold on
h2=histfit(lreturns,500,'normal');
set(h2(1),'Visible','off');
set(h2(2),'Color','g');
legend({'log-return frequencey','T locationscale'})
title('log-return relative frequency distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative frequency','fontsize',14)
set(gca,'fontsize',14)

%% complementary T cumulative distribution

pd= fitdist(lreturns(:,1),'tLocationScale');
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
legend({'pos ret','neg ret','tLocationScale'},'Location','northwest')
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
x1 = [max(lreturns1)/1000:max(lreturns1)/1000:max(lreturns1)];
loglog(x1,1-(normcdf(x1,m,s)-0.5)*2,'-m','linewidth',2)
axis([1e-5 0.5 1e-4 1])
legend({'pos ret','neg ret','normal'})
title('complemetary cumulative log-return distribution','fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
print('fat_tailsCum.eps','-depsc')

%% kstest
test=lreturns;
alpha=0.05;
pd=fitdist(test,'stable');
p1=cdf('stable',test,pd.alpha,pd.beta,pd.gam,pd.delta);
[h1,s1]=kstest(test,[test,p1],alpha)

%% ADtest
dist = makedist('stable','alpha',pd.alpha,'beta',pd.beta ,'gam',pd.gam, 'delta',pd.delta);
[h,p] = adtest(lreturns,'Distribution',dist)


%% qq plot
figure
qqplot(lreturns(:,1)),

figure,
PD = fitdist(lreturns(:,1),'GeneralizedExtremeValue'); 
qqplot(lreturns(:,1),PD)

figure,
PD = fitdist(lreturns(:,1),'tLocationScale');
qqplot(lreturns(:,1),PD)
print('qqplot d.eps','-depsc')

figure,
PD = fitdist(lreturns(:,1),'stable');
qqplot(lreturns(:,1),PD)
xlim([-2 2]);

%% fat tail power law

neg_x = sort(-lreturns(lreturns<0));
pos_x = (sort(lreturns(lreturns>0)));

[alpha1, xmin1, L1]=plfit(pos_x);%power law
[alpha2, xmin2, L1]=plfit(neg_x);

plplot(pos_x, xmin1, alpha1)
xlabel('postive log-return','fontsize',14)
ylabel('Power-law complemetary cumulative distribution','fontsize',14)
title('Compare positive log returns with power-law ','fontsize',14)
print('power1.eps','-depsc')

plplot(neg_x, xmin2, alpha2)
xlabel('negative log-return','fontsize',14)
ylabel('Power law complemetary cumulative distribution','fontsize',14)
title('Compare negative log returns with power-law ','fontsize',14)
print('power2.eps','-depsc')

%% parametric var and cvar

r_ppdf_obj = fitdist(lreturns,'stable');

%Create a list of VaRs & CVaRs where vars(i) = VaR(i percent)
par_vars = VaR(lreturns,r_ppdf_obj);
par_cvars = CVaR(lreturns,r_ppdf_obj);

%% empirical var and cvar
ys = sort(returns);
T = length(ys);

display('------------------------------')	
display('HS VaR Daily')
VaR_HS95 = -ys(int16(ceil(0.05*T)));
VaR_HS99 = -ys(int16(ceil(0.01*T)));
display(sprintf('VaR HS 95: %.9f', VaR_HS95))
display(sprintf('VaR HS 99: %.9f', VaR_HS99))
display('------------------------------')

display('------------------------------')
display('HS CVaR Daily')
ES_HS95 = -mean(ys(1:ceil(0.05*T)));
ES_HS99 = -mean(ys(1:ceil(0.01*T)));
display(sprintf('ES HS 95: %.9f', ES_HS95))
display(sprintf('ES HS 99: %.9f', ES_HS99))

ys7 = sort(returns_7);
T7 = length(ys7);

display('------------------------------')	
display('HS VaR Weekly')
VaR_HS95_7 = -ys7(int16(ceil(0.05*T7)));
VaR_HS99_7 = -ys7(int16(ceil(0.01*T7)));
display(sprintf('Weekly VaR HS 95: %.9f', VaR_HS95_7))
display(sprintf('Weekly VaR HS 99: %.9f', VaR_HS99_7))
display('------------------------------')

display('------------------------------')
display('HS CVaR Weekly')
ES_HS95_7 = -mean(ys7(1:ceil(0.05*T7)));
ES_HS99_7 = -mean(ys7(1:ceil(0.01*T7)));
display(sprintf('Weekly ES HS 95: %.9f', ES_HS95_7))
display(sprintf('Weekly ES HS 99: %.9f', ES_HS99_7))

ys14 = sort(returns_14);
T14 = length(ys14);

display('------------------------------') 
display('HS VaR Fortnightly')
VaR_HS95_14 = -ys14(int16(ceil(0.05*T14)));
VaR_HS99_14 = -ys14(int16(ceil(0.01*T14)));
display(sprintf('Fortnightly VaR HS 95: %.9f', VaR_HS95_14))
display(sprintf('Fortnightly VaR HS 99: %.9f', VaR_HS99_14))
display('------------------------------')

display('------------------------------')
display('HS CVaR Fortnightly')
ES_HS95_14 = -mean(ys14(1:ceil(0.05*T14)));
ES_HS99_14 = -mean(ys14(1:ceil(0.01*T14)));
display(sprintf('Fortnightly ES HS 95: %.9f', ES_HS95_14))
display(sprintf('Fortnightly ES HS 99: %.9f', ES_HS99_14))

ys30 = sort(returns_30);
T30 = length(ys30);

display('------------------------------')	
display('HS VaR Monthly')
VaR_HS95_30 = -ys30(int16(ceil(0.05*T30)));
VaR_HS99_30 = -ys30(int16(ceil(0.01*T30)));
display(sprintf('Monthly VaR HS 95: %.9f', VaR_HS95_30))
display(sprintf('MOnthly VaR HS 99: %.9f', VaR_HS99_30))
display('------------------------------')

display('------------------------------')
display('HS CVaR Monthly')
ES_HS95_30 = -mean(ys30(1:ceil(0.05*T30)));
ES_HS99_30 = -mean(ys30(1:ceil(0.01*T30)));
display(sprintf('Monthly ES HS 95: %.9f', ES_HS99_30))
display(sprintf('Monthly ES HS 99: %.9f', ES_HS99_30))


%% parametric var and cvar
display('------------------------------')
display('Parametric VaR')
VaR_95=portvrisk(mean(lreturns(:,1)), std(lreturns(:,1)),0.05,1);
VaR_99=portvrisk(mean(lreturns(:,1)), std(lreturns(:,1)),0.01,1);
display(sprintf('VaR 95: %.9f', VaR_95))
display(sprintf('VaR 99: %.9f', VaR_99))
display('------------------------------')

display('------------------------------')
display('Parametric CVaR')
m = mean(lreturns(:,1));
s = std(lreturns(:,1));
CVaR_95=m+s*(normpdf(-norminv(0.05,m,s),0,1))./(1-0.95);
CVaR_99=m+s*(normpdf(-norminv(0.01,m,s),0,1))./(1-0.99);
display(sprintf('CVaR 95: %.9f', CVaR_95))
display(sprintf('CVaR 99: %.9f', CVaR_99))
display('------------------------------')

%%parametric VaR weekly, fortnightly, monthly
display('------------------------------')
display('Parametric weekly VaR')
VaR_95_7=portvrisk(mean(lreturns_7(:,1)), std(lreturns_7(:,1)),0.05,1);
VaR_99_7=portvrisk(mean(lreturns_7(:,1)), std(lreturns_7(:,1)),0.01,1);
display(sprintf('VaR 95_7: %.9f', VaR_95_7))
display(sprintf('VaR 99_7: %.9f', VaR_99_7))
display('------------------------------')

display('------------------------------')
display('Parametric fortnightly VaR')
VaR_95_14=portvrisk(mean(lreturns_14(:,1)), std(lreturns_14(:,1)),0.05,1);
VaR_99_14=portvrisk(mean(lreturns_14(:,1)), std(lreturns_14(:,1)),0.01,1);
display(sprintf('VaR 95_14: %.9f', VaR_95_14))
display(sprintf('VaR 99_14: %.9f', VaR_99_14))
display('------------------------------')

display('------------------------------')
display('Parametric monthly VaR')
VaR_95_30=portvrisk(mean(lreturns_30(:,1)), std(lreturns_30(:,1)),0.05,1);
VaR_99_30=portvrisk(mean(lreturns_30(:,1)), std(lreturns_30(:,1)),0.01,1);
display(sprintf('VaR 95_30: %.9f', VaR_95_30))
display(sprintf('VaR 99_30: %.9f', VaR_99_30))
display('------------------------------')

%%parametric CVaR weekly, fortnightly, monthly
display('------------------------------')
display('Parametric weekly CVaR')
m_7= mean(lreturns_7(:,1));
s_7= std(lreturns_7(:,1));
CVaR_95_7=m_7+s_7*(normpdf(-norminv(0.05,m_7,s_7),0,1))./(1-0.95);
CVaR_99_7=m_7+s_7*(normpdf(-norminv(0.01,m_7,s_7),0,1))./(1-0.99);
display(sprintf('CVaR 95_7: %.9f', CVaR_95_7))
display(sprintf('CVaR 99_7: %.9f', CVaR_99_7))
display('------------------------------')

display('------------------------------')
display('Parametric fortnightly CVaR')
m_14 = mean(lreturns_14(:,1));
s_14 = std(lreturns_14(:,1));
CVaR_95_14=m_14+s_14*(normpdf(-norminv(0.05,m_14,s_14),0,1))./(1-0.95);
CVaR_99_14=m_14+s_14*(normpdf(-norminv(0.01,m_14,s_14),0,1))./(1-0.99);
display(sprintf('CVaR 95_14: %.9f', CVaR_95_14))
display(sprintf('CVaR 99_14: %.9f', CVaR_99_14))
display('------------------------------')

display('------------------------------')
display('Parametric monthly CVaR')
m_30 = mean(lreturns_30(:,1));
s_30 = std(lreturns_30(:,1));
CVaR_95_30=m_30+s_30*(normpdf(-norminv(0.05,m_30,s_30),0,1))./(1-0.95);
CVaR_99_30=m_30+s_30*(normpdf(-norminv(0.01,m_30,s_30),0,1))./(1-0.99);
display(sprintf('CVaR 95_30: %.9f', CVaR_95_30))
display(sprintf('CVaR 99_30: %.9f', CVaR_99_30))
display('------------------------------')

%% BOOTSTRAP VaR&CVaR 

N=200;alpha=0.01; % 1 day
a=lreturns;
T=10000;
v=zeros(T,1);
cv=zeros(T,1);
for i=1:10000
    c=randperm(numel(a));
    b=a(c(1:N));
    ys = sort(b);
    t=length(ys);
    v(i,1)=- ys(int16(ceil(alpha*N)));
    cv(i,1)=-mean(ys(1:ceil(alpha*N)));
end
var=mean(v)
cvar=mean(cv)
var_l=quantile(sort(v),0.05)
var_u=quantile(sort(v),0.95)
cvar_l=quantile(sort(cv),0.05)
cvar_u=quantile(sort(cv),0.95)
std_var_error=std(v)
std_cvar_error=std(cv)

 %% BOOTSTRAP MEAN 

randnum1=randperm(4196,1000);%产生100个随机数
randnum2=randperm(4196,1500);%

nReps = 10000;
n1 = 1500;            %sample size 1
n2 = 1000;           %sample size 2
alpha = .05;        %alpha value

x2 = lreturns(randnum1(1:1000));
x1 = lreturns(randnum2(1:1500));

myStatistic = @(x1,x2) mean(x1)-mean(x2);

sampStat = myStatistic(x1,x2);
bootstrapStat = zeros(nReps,1);
for i=1:nReps
    sampX1 = x1(ceil(rand(n1,1)*n1));
    sampX2 = x2(ceil(rand(n2,1)*n2));
    bootstrapStat(i) = myStatistic(sampX1,sampX2);
end

CI = prctile(bootstrapStat,[100*alpha/2,100*(1-alpha/2)]);

%Hypothesis test: Does the confidence interval cover zero?
H = CI(1)>0 | CI(2)<0;

clf
xx = min(bootstrapStat):.00001:max(bootstrapStat);
hist(bootstrapStat,xx);
hold on
ylim = get(gca,'YLim');
h1=plot(sampStat*[1,1],ylim,'y-','LineWidth',2);
h2=plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
h3=plot([0,0],ylim,'b-','LineWidth',2);
xlabel('Difference between means');

decision = {'Fail to reject H0','Reject H0'};
title(decision(H+1));
legend([h1,h2,h3],{'Sample mean',sprintf('%2.0f%% CI',100*alpha),'H0 mean'},'Location','NorthWest');

%% BOOTSTRAP ratio of variances

nReps = 10000;
n1 = 1000;            %sample size 1
n2 = 500;           %sample size 2

x2 = lreturns(randnum1(1:500));
x1 = lreturns(randnum2(1:1000));

myStatistic = @(x1,x2) var(x1)/var(x2);

sampStat = myStatistic(x1,x2);
CIrange = 95;
bootstrapStat = zeros(1,nReps);
for i=1:nReps
    resampX1 = x1(ceil(rand(size(x1))*length(x1)));
    resampX2 = x2(ceil(rand(size(x2))*length(x2)));
    bootstrapStat(i) = myStatistic(resampX1,resampX2);
end
CI = prctile(bootstrapStat,[50-CIrange/2,50+CIrange/2]);
disp(sprintf('Ratio of variances: %5.2f',sampStat));
disp(sprintf('%d%% Confidence interval: [%5.2f,%5.2f]',CIrange,CI(1),CI(2)));

figure
clf
xx = min(bootstrapStat):.0001:max(bootstrapStat);
hist(bootstrapStat,xx);
hold on
ylim = get(gca,'YLim');
plot(sampStat*[1,1],ylim,'y-','LineWidth',2);
plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
plot([1,1],ylim,'b-','LineWidth',2)
%set(gca,'XTick',[-10:.5:10]);
title('bootstrapping on a ratio of variances');