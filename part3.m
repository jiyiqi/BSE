clear all
load('crypto.mat')

date=crypto.Date;
price1=crypto.Close_bit;
price2=crypto.Close_eth;
price3=crypto.Close_rip;
price4=crypto.Close_nem;
price5=crypto.Close_lit;
price6=crypto.Close_das;
price7=crypto.Close_ste;
price8=crypto.Close_mon;
price9=crypto.Close_ver;
price10=crypto.Close_sia;

 price=[price1,price2,price3,price4,price5,price6,price7,price8,price9,price10];
%% logreturns and returns

returns  = diff(price)./price(1:end-1,:);
lreturns  = price2ret(price);

returns1  = diff(price1)./price1(1:end-1,:);
lreturns1  = price2ret(price1);
returns2  = diff(price2)./price2(1:end-1,:);
lreturns2  = price2ret(price2);
returns3  = diff(price3)./price3(1:end-1,:);
lreturns3  = price2ret(price3);
returns4  = diff(price4)./price4(1:end-1,:);
lreturns4  = price2ret(price4);
returns5  = diff(price5)./price5(1:end-1,:);
lreturns5  = price2ret(price5);
returns6  = diff(price6)./price6(1:end-1,:);
lreturns6  = price2ret(price6);
returns7  = diff(price7)./price7(1:end-1,:);
lreturns7  = price2ret(price7);
returns8  = diff(price8)./price8(1:end-1,:);
lreturns8  = price2ret(price8);
returns9  = diff(price9)./price9(1:end-1,:);
lreturns9  = price2ret(price9);
returns10  = diff(price10)./price10(1:end-1,:);
lreturns10  = price2ret(price10);

R =corrcoef(lreturns);
% R = corrcoef(lreturns1,lreturns2);
% polyfit(ppp,ppp1,1);
% plot(ppp,ppp1,'r+');
% 
% [b,bint,r,rint,stats]=regress(ppp,ppp1,0.05);
% b,bint,stats,rcoplot(r,rint)
% lambda=factoran(lreturns,3);
%% factor model
varname = textdata(1,2:end);% variable names
obsdates = textdata(2:end,1);
% ask 'factoran' to optimize less, otherwise there will be warnning message
optionsFactoran = statset('TolX',1e-3,'TolFun',1e-3);
% change the number of factors at which the p is larger than 0.05
[lambda,psi,T,stats,F] = factoran(lr,4,'optimopts',optionsFactoran);
stats.p;
[varname' num2cell(lambda)];
% calculate contribution of common factors to original data's variance 
Contribut = 100*sum(lambda.^2)/8; 
CumCont = cumsum(Contribut); %calculate cumulative

% when number of common factors is 3
% vbls = {'bitcoin','ethereum','ripple','nem','litecoin','stellar','dash','monero','verge','siacon'};
% biplot(lambda,'LineWidth',2,'MarkerSize',20,'varlabels',vbls)