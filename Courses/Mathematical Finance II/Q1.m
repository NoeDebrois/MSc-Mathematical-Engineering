%% Q.1.1 Calibration of a local volatility model with ECORP market data :
% Reset
clc;
clearvars;
%% Q.1.1 Load market data from MarketData.xlsx
filePath = '/Users/noedebrois/Desktop/Politecnico/Mathematical Finance/Local Volatility Model/MarketData.xlsx';

%% Q.1.1 Extract E-CORP market data :
% EXTRACT MARKET EXPIRIES :
cellRangeT = 'B6:I6'; 
T = xlsread(filePath, cellRangeT);

% EXTRACT FORWARD AT MARKET EXPIRIES :
cellRangeFWD = 'B8:I8'; 
Fwd = xlsread(filePath, cellRangeFWD);

% EXTRACT MARKET STRIKES :
cellRangeK = 'B24:I30'; 
K = xlsread(filePath, cellRangeK);

% EXTRACT MARKET IMPLIED VOLATILITY :
cellRangeMKTVOL = 'B13:I19'; 
MktVol = xlsread(filePath, cellRangeMKTVOL);

%% Q.1.1 Calibration :
% Normalize market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Solver settings
Lt = 10;
Lh = 200;
K_min = 0.1;
K_max = 3;
Scheme = 'cn';

% Calibration settings
Threshold = 0.0010;
MaxIter = 100;

% Call the calibration function
[V, ModelVol, MaxErr] = calibrator(T, K_norm, MktVol, Threshold, MaxIter, Lt, Lh, K_min, K_max, Scheme);
% V is the calibrated local volatility matrix
% ModelVol is the model implied volatilities at the normalized strikes
% MaxErr is the vector of calibration error at each iteration

%% Q.1.1 Plots
% Plot local volatility function vs market implied volatility
figure;
plot(K(:, 1), MktVol(:, 1), 'o', K(:, 1), ModelVol(:, 1), ':.', K(:, 1), V(:, 1), ':.b', 'linewidth', 1.5);
title('Calibrated model and local volatility for asset E CORP');
legend('MktVol', 'ModelVol', 'LocalVol');

figure;
plot(MaxErr, '.', 'MarkerSize', 15);
title('Calibration error at each iteration of the fixed-point calibration');

%% Q.1.2 Pricing of two call options with the previous model & Dupire or MC :
% EXTRACT DISCOUNT FACTOR AT MARKET EXPIRIES
cellRangeDF = 'B7:I7'; 
DF = xlsread(filePath, cellRangeDF);

% EXTRACT SPOT PRICE :
Spot = xlsread(filePath, 'B4:B4');

%% Q.1.2 Calibration of r and q :
[r,q] = calibrate_r_q(Spot,T,DF,Fwd);

%% Q.1.2 Monte Carlo pricing of 2 call options using the above model :
% Option data
expiry = 0.5; % see Q.1.2
Kappa_val = [0.9 1.1]; % see Q.1.2

N = 1000000; % MC simulations
M = 50; % timesteps

% MC simulation and option pricing
for i = 1:length(Kappa_val)
    Kappa = Kappa_val(i);
    strike = Kappa * Spot;

    % MC simulation (LV)
    S = lv_simulation_log(T, Spot, r, q, V, K, N, M, expiry);

    % Option price (LV)
    discount_factor = discount(T, r, expiry);
    lv_price = discount_factor * mean(max(S(1, :) - strike, 0));
    
    % LV implied spot volatility
    fwd = forward(Spot,T,r,q,expiry);
    impl_vol_lv = blsimpv(fwd,strike,0,expiry,lv_price/discount_factor);
    
    % Display the result
    fprintf('Monte Carlo LV Price for Kappa = %.2f: %.4f\n', Kappa, lv_price);
    fprintf('Monte Carlo LV Implied Volatility for Kappa = %.2f: %.4f\n', Kappa, impl_vol_lv);
end

%% Q.1.3 Monte Carlo pricing of 2 fwd starting options :
% Forward starting option data
start_date = 2;
expiry_date = 2.5;

% MC simulation for forward starting options
for i = 1:length(Kappa_val)
    Kappa = Kappa_val(i);

    % Calculate start and expiry times
    T1 = start_date;
    T2 = expiry_date;
    
    % Generate paths using LV simulation
    S1 = lv_simulation_log(T, Spot, r, q, V, K, N, M, T1);
    S2 = lv_simulation_log(T, Spot, r, q, V, K, N, M, T2);

    % Calculate the option price
    discount_factor_T1 = discount(T, r, T1);
    discount_factor_T2 = discount(T, r, T2);
    P = discount_factor_T2 * mean(max(S2 - Kappa * S1, 0));

    % Compute implied forward volatility
    fwd_start = forward(Spot, T, r, q, T1);
    fwd_end = forward(Spot, T, r, q, T2);
    model_impl_fwd_volatility = blsimpv(fwd_end, Kappa * fwd_start, 0, T2 - T1, P / discount_factor_T2);

    % Display the result
    fprintf('\nForward Starting Option MC price for Kappa = %.2f:  %.4f\n', Kappa, P);
    fprintf('MC Implied Forward Volatility: %.4f', model_impl_fwd_volatility);
end

%% Local volatility surface :
qt = 2:.05:3;
qk = 4000:500:10000;
s=zeros(length(qt),length(qk));
for a=1:length(qt)
    for c=1:length(qk)
        s(a,c) = localvol(T,K,V,qt(a),qk(c));
    end
end

figure;
surf(s);
xlabel('Strikes');
ylabel('Time to Expiry');
zlabel('Local Volatility');
title({'Local Volatility surface:';'piecewise-constant interpolation on time';'linear interpolation on strikes with flat extrapolation'});

%% Q.1.4 Spot smile and forward smile :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING SPOT START CALL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% spot start option data
Expiry = 0.5;

% model forwards at T
Fwd = zeros(1,length(T));
for i=1:length(T)
   Fwd(i) = forward(Spot,T,r,q,T(i));
end

% normalized LV strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Dupire solver details
Lt = 10;
Lh = 200;
K_min = 0.1;
K_max = 3;
Scheme = 'cn';

% compute price of spot-start call options 
[ k, C ] = solve_dupire( T, K_norm, V, Expiry, Lt, Lh, K_min, K_max, Scheme);    

% compute model implied volatilities
perc_strikes_spot_start=[];
model_impl_vol_spot_start=[];
for i=1:length(k)
    if k(i)>0.3 && k(i)<1.5
        perc_strikes_spot_start = [perc_strikes_spot_start k(i)];
        model_impl_vol_spot_start = [model_impl_vol_spot_start blsimpv(1,k(i),0,Expiry,C(i))];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING FORWARD START CALL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% forward-start option data
expiry(1) = 2;
expiry(2) = expiry(1) + 0.5;

% additional market data
discount_factor = discount(T,r,expiry(2));
fwd(1) = forward(Spot,T,r,q,expiry(1));
fwd(2) = forward(Spot,T,r,q,expiry(2));

N=1000000; %MC simulations
M=100; %timesteps

% MC simulation
S = lv_simulation_log(T,Spot,r,q,V,K,N,M,expiry);

% option prices
perc_strikes = 0.3:0.05:1.5;
model_impl_vol=[];
for x = perc_strikes
    P = discount_factor*mean(max(S(2,:) - x*S(1,:),0));
    model_impl_vol = [model_impl_vol, blsimpv(fwd(2),x*fwd(1),0,expiry(2)-expiry(1),P/discount_factor)];
end

plot(perc_strikes_spot_start,model_impl_vol_spot_start,perc_strikes,model_impl_vol);
title('Model implied vol with option maturity 6m');
legend('Spot impl vol','Fwd impl vol (starts in 2y)')

%% Q.1.4 Compute skew of spot smile and forward smile :
% Calculate skewness
skewness_spot_smile = skewness(model_impl_vol_spot_start);
skewness_fwd_smile = skewness(model_impl_vol);

% Display the result
fprintf('\nSkewness of the spot smile: %.4f\n', skewness_spot_smile);
fprintf('Skewness of the fwd smile: %.4f\n', skewness_fwd_smile);





















