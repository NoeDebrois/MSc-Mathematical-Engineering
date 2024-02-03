%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap the Zero-Coupon curve from a 'strip' of Interest Rate Swaps   %
% See Bjork 4th Ed. Ch. 19                                                %
% AN 10/11/2023                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plots visible also to students sitting in the last row
set(0,'defaultAxesFontSize',16)

%% Pricing functions
p = @(y,t,T) exp(-y.*(T-t));             % Zero coupon price
y = @(p,t,T) - log(p)./(T-t);            % Zero coupon yield

%% Market data: quotes of the strip of IRS (BID) 
%  Pay LIBOR quarterly - receive FIX semiannual (bond basis: 30/360)
IRS =    [   1.00 5.677; 2.00 5.265; 3.00 5.019; 4.00 4.867; 5.00 4.885; 
             6.00 4.865; 7.00 4.8415; 8.00 4.824; 9.00 4.8427;10.0 4.8463;
             11.0 4.849; 12.0 4.8542; 13.0 4.861; 14.0 4.855; 15.0 4.851;
             20.0 4.797; 25.0 4.709; 30.0 4.6215; 40.0 4.412; 50.0 4.202];
freq = 2;               % Semiannual fixed coupons
t = 0;                  % Issue date of the par-yield bonds (now)

%% The first node of the ZC curve is the 'old good yield-to-maturity'
i = 1;
Tmat = IRS(i,1);    % Maturity of the i-th par-yield bond
cpn = IRS(i,2);     % Coupon rate of the i-th par-yield bond
guess = cpn;
fun = @(y) price([t y; Tmat y],t,Tmat,cpn,freq)-1.0;
YTM = fzero(fun,guess);      %Yield to Maturity (percentage points)
zCurve = [Tmat YTM]; % Bootstrap baseline: the first node of the ZC curve

%% Loop across the IRS strip and bootstrap the curve
for i = 2:height(IRS)   % Lenght of the IRS strip -> nodes of the ZC curve
    Tmat = IRS(i,1);    % Maturity of the i-th par-yield bond
    cpn = IRS(i,2);     % Coupon rate of the i-th par-yield bond
    guess = cpn;
    fun = @(y) price([zCurve; Tmat y],t,Tmat,cpn,freq)-1.0;
    yi = fzero(fun,guess);      
    zCurve = [zCurve; Tmat yi];
end

%% Term structure of the Zero-Coupon bond prices
ZC_nodes  = zCurve(:,1)';
ZC_rates  = zCurve(:,2)'/100;
ZC_prices = p(ZC_rates,t,ZC_nodes);
pCurve = [ZC_nodes;ZC_prices]'; 

%% Instantaneous forwards
fwd_rates = -df(ZC_nodes,log(ZC_prices))*100;
fCurve = [ZC_nodes;fwd_rates]';        % The instantaneous forward curve

%% Plot curves
figure
plot(zCurve(:,1),zCurve(:,2),'x','Linewidth', 2,'LineStyle','-')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Zero Coupon Curve')
xlabel('Time (year frac)')

figure
plot(pCurve(:,1),pCurve(:,2),'x','Linewidth', 2,'LineStyle','-')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Zero Coupon bond prices')
xlabel('Time (year frac)')

figure
plot(fCurve(:,1),fCurve(:,2),'x','Linewidth', 2,'LineStyle','-')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Instantaneous Forward Curve')
xlabel('Time (year frac)')

%% Test
%Ten years par-yield bond
calculationDate = '25 oct 2023';
maturityDate= '25 oct 2033';
Tmat = yearfrac(calculationDate, maturityDate,5);
P_mkt = 100.00;                            
P_model = 100*price(zCurve,t,Tmat,4.8463,2) %Price from discounted cash-flows