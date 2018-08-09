function [ CVaRs ] = CVaR(data_set, pd_obj)
%Calculate expected shortfalls of data set from given distribution object
    %Note that CVaRs(i) is loss of i% or more, NOT confidence level i%

%Error = acceptable error either side of target percent (/100)
error = 0.005;

%Find edges of cdf for data from the distribution object
x = linspace(-abs(min(data_set)),abs(max(data_set)),2);
data_cdf = cdf(pd_obj,x);

%Make sure cdf has values from error/100 to 1-error/100
low_m = 1;
high_m = 1;
while min(data_cdf) > error/100
    low_m = low_m+1;
    x = linspace(-low_m*abs(min(data_set)),abs(max(data_set)),2);
    data_cdf = cdf(pd_obj,x);
end
while max(data_cdf) < (1-error/100)
    high_m = high_m+1;
    x = linspace(-low_m*abs(min(data_set)),high_m*abs(max(data_set)),2);
    data_cdf = cdf(pd_obj,x);
end

%Set initial guess of good cdf
x = linspace(-low_m*abs(min(data_set)),high_m*abs(max(data_set)),1000);
data_cdf = cdf(pd_obj,x);
    
%Make sure will find values in cdf within error of any given percent
if max(diff(data_cdf)) > 2*error
    x = linspace(-low_m*abs(min(data_set)),high_m*abs(max(data_set)),10000);
    data_cdf = cdf(pd_obj,x);
end
if max(diff(data_cdf)) > 2*error
    x = linspace(-low_m*abs(min(data_set)),high_m*abs(max(data_set)),100000);
    data_cdf = cdf(pd_obj,x);
end

%Trim cdf to more manageable size based on acceptable error
if sum(data_cdf<error/100) ~= 0
    x = x(sum(data_cdf<error/100):end);
    data_cdf = data_cdf(sum(data_cdf<error/100):end);
end
if sum(data_cdf>1-error/100) ~=0
    x = x(1:end-sum(data_cdf>1-error/100)+1);
    data_cdf = data_cdf(1:end-sum(data_cdf>1-error/100)+1);
end

%Calculate VaR from cdf
VaR_calc = @(percent) mean(x(percent/100-error<data_cdf & ...
    data_cdf<percent/100+error));

CVaRs = zeros(1,100);
vars = zeros(1,1001);
for i = 1:100
    %Calculate VaR from VaR(0) to VaR(percent) in 1001 increments
    for j = linspace(0,i,1001)
        index = round(j*1000/i+1);
        vars(index) = VaR_calc(j);
    end

    %Integrate over computed VaRs to find average
    CVaRs(i) = trapz(vars)/1000; %trapz b/c integration
    
end

end

