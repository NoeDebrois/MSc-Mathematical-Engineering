%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yield-to-Maturity and Duration of a coupon bond                         %
% See Bjork 4th Ed. Ch. 19                                                %
% AN 13/10/2023                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pricing functions
p = @(y,t,T) exp(-y.*(T-t));             % Zero coupon price
price = @(y,t,T,CF) sum(p(y,t,T) .* CF); % Coupon bond price

%% Define the timeframe
%Ten years (original maturity) T-bond
issueDate = '15 aug 23';
calculationDate = '13 oct 2023';
t = yearfrac(issueDate, calculationDate)
%% Pricing parameters
delta = 0.5;        %Semiannual coupons
K = 100;            %Principal
T1 = delta;         %Next coupon date (yearfrac)
Tn = 10;            %Last coupon date (yearfrac)
n = Tn/delta;       %Number of coupons
r = 3.875/100;      %Coupon Rate
%% Accrued interest
accruedInterest = K*r*t;
%%  Coupon's schedule and cash-flows (see Eq. 19.16 and Eq. 19.17)
T = linspace(T1,Tn,n);
CF = K*ones(1,n)*r*delta;
CF(1,n) = CF(1,n) + K;
%% Market data
P = 93.164;                 %Clean price from WSJ page
P = P +accruedInterest;     %Dirty price
%% Yield to Maturity (percentage points)
guess = 0.05;
fun = @(y) price(y,t,T,CF)-P;
YTM = fzero(fun,guess)*100
%% Duration
D = sum(T.*CF.*p(YTM/100,t,T))/price(YTM/100,t,T,CF)