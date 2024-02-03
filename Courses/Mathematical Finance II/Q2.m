%% Q.2
% Reset
clc;
clearvars;
%% Load market data from MarketData XY.xls
filePath = '/Users/noedebrois/Desktop/Politecnico/Mathematical Finance/Local Volatility Model/MarketData.xlsx';
sheetName = 'FAIL';
%% Extract data
% EXTRACT MARKET EXPIRIES :
cellRangeT = 'B6:C6'; 
T = xlsread(filePath, sheetName, cellRangeT);

% EXTRACT FORWARD AT MARKET EXPIRIES :
cellRangeFWD = 'B8:C8'; 
Fwd = xlsread(filePath, sheetName, cellRangeFWD);

% EXTRACT MARKET STRIKES :
cellRangeK = 'B22:C26'; 
K = xlsread(filePath, sheetName, cellRangeK);

% EXTRACT MARKET IMPLIED VOLATILITY :
cellRangeMKTVOL = 'B13:C17'; 
MktVol = xlsread(filePath, sheetName, cellRangeMKTVOL);

% EXTRACT SPOT PRICE :
Spot = xlsread(filePath, sheetName, 'B4:B4');

% EXTRACT DISCOUNT FACTOR :
DiscFact = xlsread(filePath, sheetName, 'B7:C7');
%% Computations : THEY FAIL !
% % Normalize market strikes
% K_norm = K ./ Fwd;
% 
% % Dupire solver settings
% Lt = 10;
% Lh = 200;
% K_min = 0.1;
% K_max = 3;
% Scheme = 'cn';
% 
% % Calibration settings
% Threshold = 0.0010;
% MaxIter = 100;
% 
% % Call the calibration function
% [V, ModelVol, MaxErr] = calibrator(T, K_norm, MktVol, Threshold, MaxIter, Lt, Lh, K_min, K_max, Scheme);
% 
% % Plot local volatility function vs market implied volatility
% figure;
% plot(K(:, 1), MktVol(:, 1), 'o', K(:, 1), ModelVol(:, 1), ':.', K(:, 1), V(:, 1), ':.b', 'linewidth', 1.5);
% title('Calibrated model and local volatility for asset E CORP');
% legend('MktVol', 'ModelVol', 'LocalVol');
% 
% figure;
% plot(MaxErr, '.', 'MarkerSize', 15);
% title('Calibration error at each iteration of the fixed-point calibration');
%% Let's check that option prices are a convex function of the strike prices :
% Calibration of r and q from the model :
[r,q] = calibrate_r_q(Spot,T,DiscFact,Fwd);
%%
% Pricing : 
T1 = T(:,1);
K1 = K(:,1);
r1 = r(:,1);
q1 = q(:,1);
MktVol1 = MktVol(:,1);
priceT1 = blsprice(Spot,K1,r1,T1,MktVol1,q1);

T2 = T(:,2);
K2 = K(:,2);
r2 = r(:,2);
q2 = q(:,2);
MktVol2 = MktVol(:,2);
priceT2 = blsprice(Spot,K2,r2,T2,MktVol2,q2);
%% Plot :
figure;
plot(K1,priceT1);
hold;
plot(K2,priceT2);
title('Black-Scholes price as function of volatility for fixed K,T');
