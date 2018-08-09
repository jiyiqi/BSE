clear
close all

%%%%%%%%%%%%%%%
% load high frequency prices for Apple 
load('AAPL_1sec_3Dic2012.mat') % data from http://www.dukascopy.com/swiss/english/marketwatch/historical/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find range with active trading
t1=find(diff(ClosingPrice)~=0,1,'first')+1;
t2=find(diff(ClosingPrice)~=0,1,'last');
Dates = Dates(t1:t2);
pr=ClosingPrice(t1:t2); 
%%%%%%%%%%%%%%%%%%%
% stylized fact -6-(seconds)
%%%%%%%%%%%%%%%%%%%
syb = 'v<>^s';
col = 'rbmgc'
k=0;
mm=[1 10 30 60 300];
figure
for m=mm
    k=k+1;
    lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
    [f,b]=hist(lr,20);
    semilogy(b,f/max(f),['-',syb(k),col(k)])
    hold on
end
axis([-0.002 +0.002 1e-4 1.1])
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('lor-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('AAPL 1 sec data 3 Dic 2012')
print('1secRetDistr.eps','-depsc')
%
h=GenHurst(pr,2);
k=0;
figure
for m=mm
    k=k+1;
    lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
    [f,b]=hist(lr,20);
    semilogy(b/m.^h,f/max(f),['-',syb(k),col(k)])
    hold on
end
text(0.0005,0.05,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.0015 +0.0015 1e-4 1.1])
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('AAPL 1 sec data 3 Dic 2012')
print('1secRetDistrScaled.eps','-depsc')


returns = diff(log(pr));
%%%%%%%%%%%%%%%%%%%
% stylized fact -2-
%%%%%%%%%%%%%%%%%%%
figure
[a1,lags1] = autocorr(returns,150);
plot(lags1,a1,'-r')
hold on
[a2,lags2] = autocorr(abs(returns),150);
plot(lags2,a2,'-m')
[a3,lags3] = autocorr(returns.^2,150);
plot(lags3,a3,'-b')
plot([0 300],[0 0],'-k')
axis([0 150 -.25 0.25])
xlabel('lags (seconds)','fontsize',14)
ylabel('autocorrelation','fontsize',14)
title('AAPL autocorrelation of returns','fontsize',14)
set(gca,'fontsize',14)
legend({'returns','|returns|','returns^2'})
print('autocorrelation1sec.eps','-depsc')


