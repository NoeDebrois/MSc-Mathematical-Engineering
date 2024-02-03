function [price] = price(ZCcurve,t,Tmat,cpn,freq)
%   Price of the coupon bond by discounted cashflows
%   ZC_curve: Table of ZC rates (cont. comp. 30/360)
%       Column #1: Tj (year frac)
%       Column #2: y(t,Tj) (percentage points)
%   t:        Calculation time
%   Tmat:     Bond maturity (year frac)
%   cpn:      Coupon rate (percentage points)
%   freq:     Payment frequency
%        1: annual
%        2: semi-annual
%        4: quarterly
%       12: monthly

%% Pricing function
p = @(y,t,T) exp(-y.*(T-t));             % Zero coupon price

%% The interest rate curve is allocated into two vectors
ZC_nodes = ZCcurve(:,1);
ZC_rates = ZCcurve(:,2)/100;

%% Pricing parameters
K = 1;                  %Principal
delta = 1/freq;         %Time intervals
n = ceil(Tmat/delta);   %Number of coupons
Tn = Tmat;              %Last coupon date (yearfrac)
T1 = Tn-delta*(n-1);    %Next coupon date (yearfrac)
r_cpn = cpn/100;        %Coupon Rate

%%  Coupon's schedule and cash-flows (see Eq. 19.16 and Eq. 19.17)
T = linspace(T1,Tn,n);
CF = K*ones(1,n)*r_cpn*delta;
CF(1,n) = CF(1,n) + K;

%% Accrued interest assuming t=0 is the last coupon date (or the issue date) of the bond
accruedInterest = K*r_cpn*t;

%% ZC rates at the payment dates y(t,Ti)
r = ZC_rates(1);                    %Approximated Short-term rate
y = interp1(ZC_nodes,ZC_rates,T,'linear',r);
% y = interp1(ZC_nodes,ZC_rates,T,'spline',r);

%% Coupon bond price
price = sum(p(y,t,T) .* CF)+ accruedInterest;

end