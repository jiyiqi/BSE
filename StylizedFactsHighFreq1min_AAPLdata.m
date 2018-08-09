close all
clear

%%%%%%%%%%%%%%%
% load high frequency prices for Apple 
load('AAPL_1min_Nov2012.mat') % data from http://www.dukascopy.com/swiss/english/marketwatch/historical/
[~,ui] = unique(Day);
ud = Day(sort(ui));
cumVols = zeros(length(Volume(Day==ud(2))),1);
mm = [1 2 3 4 5 6 7 8 9 10 12 15 20 25 30 35 40 45 50 55 60]
rets= [];
vol = [];
Vol = [];
for m=mm,lr{m} = [];end
for d=ud(2:end)' %nb first day is shorter
    dd = find(Day==d );
    ddd = find(Day==d & Volume > 1e5);
    hh1= find(Hour(ddd)==15);
    mm1= find(Minute(hh1)==50);
    t1 = ddd(hh1(mm1));
    hh2= find(Hour(ddd)==21);
    mm2= find(Minute(hh2)==10);
    t2 = ddd(hh2(mm2));
    pr = ClosingPrice(t1:t2);
%     plot(pr)
%     pause(1)
%     clf
    for m=mm
        lr{m} = [lr{m};log(pr((m+1):m:end))-log(pr(1:m:(end-m)))];
    end
    rets = [rets;diff(pr(:))];
    m = 150;
    x  = diff(pr(:));
    Y1 = cumsum(x,1) ;
    Y2 = cumsum(x.^2,1) ;
    volatility = sqrt((Y2((m+1):end)-Y2(1:(end-m)))/m-((Y1((m+1):end)-Y1(1:(end-m)))/m).^2) ;
    vol = [vol;volatility(:)];
    x  = Volume(t1:(t2-1));
    Y1 = cumsum(x,1) ;
    Y2 = cumsum(x.^2,1) ;
    cVol = sqrt((Y2((m+1):end)-Y2(1:(end-m)))/m-((Y1((m+1):end)-Y1(1:(end-m)))/m).^2) ;
    Vol = [Vol;cVol(:)];
    cumVols  = cumVols + Volume(dd);
end
pri = cumsum(rets); % an artificial trading only price
figure
plot(pri)
h=GenHurst(pri,2);

%%%%%%%%%%%%%%%%%%%
% stylized fact -3-(minutes)
c=corrcoef(Vol,vol);
figure
plot(Vol,vol,'+')
text(max(Vol)*.1,max(vol)*.9,['correlation coeff =',num2str(c(1,2))],'fontsize',14)
xlabel('volume','fontsize',14)
ylabel('volatility','fontsize',14)
set(gca,'fontsize',14)
print('volume_volatility1min.eps','-depsc')

figure
subplot(2,1,1)
plot(vol)
xlabel('date','fontsize',14)
ylabel('volatility','fontsize',14)
set(gca,'fontsize',14)
subplot(2,1,2)
plot(Vol)
xlabel('date','fontsize',14)
ylabel('volume','fontsize',14)
set(gca,'fontsize',14)
print('volatility1min.eps','-depsc')

%%%%%%%%%%%%%%%%%%%
% stylized fact -4-(minutes)
%%%%%%%%%%%%%%%%%%%
figure
plot(Dates(dd),[cumVols])
datetick
set(gca,'fontsize',14)
xlabel('time','fontsize',14)
ylabel('volume','fontsize',14)
title('AAPL  cumulate intraday volumes Nov 2012')
print('seasonality1min.eps','-depsc')

%%%%%%%%%%%%%%%%%%%
% stylized fact -6-(minutes)
%%%%%%%%%%%%%%%%%%%
h=GenHurst(pri,2);
bins=[];
bins1=[];
freqs=[];
for m=mm
    s(m) = std(lr{m});
    [f,b]=hist(lr{m},20);
    bins = [bins,b'/m];
    bins1 = [bins1,b'/m^h];
    freqs = [freqs,f'/max(f)];
end
figure
semilogy(bins,freqs,'-+')
%text(0.002,0.2,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.006 +0.006 0 1.1])
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('AAPL 1 min data Nov 2012')
print('RetDistr1min.eps','-depsc')

figure
semilogy(bins1,freqs,'-+')
hold on
x=-0.01:0.0001:0.01;
y=normpdf(x,0,s(1));
semilogy(x,y/max(y)*0.55,'--k')
text(0.002,0.2,['H(2) = ',num2str(h)],'fontsize',14)
axis([-0.006 +0.006 5e-4 1.1])
legend(num2str(mm'))
set(gca,'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('relative freq','fontsize',14)
title('AAPL 1 min data Nov 2012')
print('RetDistrScaled1min.eps','-depsc')

