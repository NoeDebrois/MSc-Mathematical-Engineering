%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Zero-coupon curve                                                   %
% See Bjork 4th Ed. Ch. 19                                                %
% AN 27/10/2023                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pricing functions
p = @(y,t,T) exp(-y.*(T-t));             % Zero coupon price
y = @(p,t,T) - log(p)./(T-t);            % Zero coupon yield

%% Test data
pMkt =    [ 0.25 0.9880; 0.50 0.9748; 0.75 0.9612; 1.00 0.9471; 
             2.00 0.9048; 3.00 0.8711; 4.00 0.8403; 5.00 0.8106; 
             6.00 0.7815; 7.00 0.7542; 8.00 0.7273; 9.00 0.7008; 
             10.0 0.60];

t = 0;
x = pMkt(:,1)';
y = (y(pMkt(:,2),t,pMkt(:,1))*100)';
zCurve = [x;y]';            % The zero-coupon curve

%%
f = -df(x,log(pMkt(:,2))')*100;
fCurve = [x;f]';            % The instantaneous forward curve

%%
% Plot curves
figure
plot(zCurve(:,1),zCurve(:,2),'x','Linewidth', 2,'LineStyle','-')
%grid
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Zero Coupon Curve')
xlabel('Time (year frac)')

%%
figure
plot(fCurve(:,1),fCurve(:,2),'x','Linewidth', 2,'LineStyle','-')
%grid
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Instantaneous Forward Curve')
xlabel('Time (year frac)')

%% Define the timeframe
%Ten years (original maturity) T-bond
calculationDate = '13 oct 2023';
maturityDate= '15 aug 2033';

Tmat = yearfrac(calculationDate, maturityDate);

%% test
P_mkt = 93.164;                             %Clean price from WSJ page
P_model = 100*price(zCurve,t,Tmat,3.875,2) %Clean price from discounted cash-flows