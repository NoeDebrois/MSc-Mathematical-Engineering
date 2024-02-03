%% Q.3.1 to Q.3.4
%% Q.3.1 Market data of the FX asset Y = EUR/USD :
% Reset :
clc;
clearvars;
% Paths :
filePath = '/Users/noedebrois/Desktop/Politecnico/Mathematical Finance/Local Volatility Model/MarketData.xlsx';
sheetName = 'EURUSD';
% EXTRACT SPOT PRICE :
Spot = xlsread(filePath, sheetName, 'B4:B4');
% EXTRACT MARKET EXPIRIES (year fractions) :
cellRangeT = 'B6:H6'; 
T = xlsread(filePath, sheetName, cellRangeT);
% EXTRACT DOMESTIC DISCOUNT FACTOR AT MARKET EXPIRIES
cellRangeDomDF = 'B7:H7'; 
DomDF = xlsread(filePath, sheetName, cellRangeDomDF);
% EXTRACT FOREIGN DISCOUNT FACTOR AT MARKET EXPIRIES
cellRangeForDF = 'B8:H8'; 
ForDF = xlsread(filePath, sheetName, cellRangeForDF);
% EXTRACT FORWARD AT MARKET EXPIRIES :
cellRangeFWD = 'B9:H9'; 
Fwd = xlsread(filePath, sheetName, cellRangeFWD);
% EXTRACT MARKET IMPLIED VOLATILITY :
cellRangeMKTVOL = 'B14:H18'; 
MktVol = xlsread(filePath, sheetName, cellRangeMKTVOL);
% EXTRACT DELTAS :
cellRangeDeltas = 'A14:A18';
Deltas = xlsread(filePath, sheetName, cellRangeDeltas);

%% Q.3.1 Computation of the market strikes {K_i,j} from the quoted deltas :
K = zeros(length(Deltas), length(T));
for i = 1:length(Deltas)
    for j = 1:length(T)
        K(i,j) = fzero(@(Strike) blsdelta(Fwd(j),Strike,0,T(j),MktVol(i,j))-(1-Deltas(i)/100), Fwd(j));
    end
end

%% Q.3.1 Calibration of the local volatility model (max_calib_err = 0.001) :
% normalized market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Dupire solver settings
Lt = 20;
Lh = 300;
K_min = 0.5;
K_max = 2.5;
Scheme = 'cn';

% calibration settings
Threshold = 0.0010;
MaxIter = 100;

[V, ModelVol, MaxErr] = calibrator(T,K_norm,MktVol,Threshold,MaxIter,Lt,Lh,K_min,K_max,Scheme);

% plot local volatility function vs market implied volatility
figure;
plot(K(:,1),MktVol(:,1),'o',K(:,1),ModelVol(:,1),':.',K(:,1),V(:,1),':.b','linewidth',1.5);
title('Calibrated model and local volatility for asset EUR/USD');
legend('MktVol','ModelVol','LocalVol');
xlabel('\{K_{i,1}\}_i \in \{1, 2, 3, 4, 5\}');
ylabel('Mkt, Model and Local volatilities');

figure;
plot(MaxErr,'.','MarkerSize',15);
title('calibration error at each iteration of the fixed-point calibration');

%% Q.3.2 Monte-Carlo prices for 1 plain vanilla & 1 digital options :
N = 100000;
M = 100;

% Calibration or r and q (with Foreign Discount Factor ! ??) :
[rf,qf] = calibrate_r_q(Spot,T,ForDF,Fwd);

% Option data
T5 = xlsread(filePath, sheetName, 'F6:F6');
K5_2 = K(2,5); % The strike computed at 3.1, corresponding to 25-Delta

% MC simulation (Local Vol)
S = lv_simulation_log(T,Spot,rf,qf,V,K,N,M,T5); % scenari

% Plain Vanilla Option with T_5 and K_{5,2} :
discount_factor = discount(T,rf,T5);
% LV price of a call option
price_lv_vanilla_vec = discount_factor*max(S(1,:) - K5_2,0);
vanilla_option_price = mean(price_lv_vanilla_vec); % valuation
% price computed using (rf,qf) (is it good?)

% Digital Option with T_5 and K_{5,2} :
% LV price of a digital option
price_lv_digital_vec = discount_factor*(S(1,:) > K5_2); % payoff is 1 if S_T > K_{5,2}
digital_option_price = mean(price_lv_digital_vec); % valuation
% again, price computed using (rf,qf) (is it good?)

%% Q.3.2 Confidence interval of level 95% for the option prices :
confidence_level = 0.95;  % Desired confidence level

% Compute Sample Mean and Sample Variance
mean_plain_vanilla = vanilla_option_price;
variance_plain_vanilla = std(price_lv_vanilla_vec)^2;

mean_digital_option = digital_option_price;
variance_digital_option = std(price_lv_digital_vec)^2;

% Calculate Monte Carlo Error
monte_carlo_error_plain_vanilla = sqrt(variance_plain_vanilla / N);
monte_carlo_error_digital_option = sqrt(variance_digital_option / N);

z_score = norminv((1 + confidence_level) / 2, 0, 1);

% Compute Confidence Interval for Plain Vanilla Option
lower_bound_plain_vanilla = mean_plain_vanilla - z_score * monte_carlo_error_plain_vanilla;
upper_bound_plain_vanilla = mean_plain_vanilla + z_score * monte_carlo_error_plain_vanilla;

% Compute Confidence Interval for Digital Option
lower_bound_digital_option = mean_digital_option - z_score * monte_carlo_error_digital_option;
upper_bound_digital_option = mean_digital_option + z_score * monte_carlo_error_digital_option;

% Display Results
fprintf('### MONTE CARLO - LOCAL VOLATILITY MODEL ###\n')
fprintf('Plain Vanilla Option:\n');
fprintf('Monte Carlo Price: %f\n', vanilla_option_price);
fprintf('Confidence Interval (%.2f%%): [%.6f, %.6f]\n', confidence_level * 100, lower_bound_plain_vanilla, upper_bound_plain_vanilla);

fprintf('\nDigital Option:\n');
fprintf('Monte Carlo Price: %f\n', digital_option_price);
fprintf('Confidence Interval (%.2f%%): [%.6f, %.6f]\n', confidence_level * 100, lower_bound_digital_option, upper_bound_digital_option);

%% Q.3.3 Monte-Carlo prices and confidence intervals with Black dynamics
sigma = MktVol(2,5); % 25-Delta volatility at expiry T5

% MC simulation (Black dynamics)
S_B = black_simulation_log(T,Spot,rf,qf,sigma,N,T5); % scenari

% MC price of a plain vanilla option with T_5 and K_{5,2} under Black
% dynamics :
price_B_vanilla_vec = discount_factor*max(S_B(1,:) - K5_2,0);
vanilla_option_B_price = mean(price_B_vanilla_vec); % valuation

% MC price of a digital option with T_5 and K_{5,2} under Black dynamics :
price_B_digital_vec = discount_factor*(S_B(1,:) > K5_2);
digital_option_B_price = mean(price_B_digital_vec); % valuation
% again, price computed using (rf,qf) (is it good?)

%% Q.3.3 Confidence interval of level 95% for the option prices under Black dynamics :
% Compute Sample Mean and Sample Variance
mean_plain_vanilla_B = vanilla_option_B_price;
variance_plain_vanilla_B = std(price_B_vanilla_vec)^2;

mean_digital_option_B = digital_option_B_price;
variance_digital_option_B = std(price_B_digital_vec)^2;

% Calculate Monte Carlo Error
monte_carlo_error_plain_vanilla_B = sqrt(variance_plain_vanilla_B / N);
monte_carlo_error_digital_option_B = sqrt(variance_digital_option_B / N);

% Compute Confidence Interval for Plain Vanilla Option
lower_bound_plain_vanilla_B = mean_plain_vanilla_B - z_score * monte_carlo_error_plain_vanilla_B;
upper_bound_plain_vanilla_B = mean_plain_vanilla_B + z_score * monte_carlo_error_plain_vanilla_B;

% Compute Confidence Interval for Digital Option
lower_bound_digital_option_B = mean_digital_option_B - z_score * monte_carlo_error_digital_option_B;
upper_bound_digital_option_B = mean_digital_option_B + z_score * monte_carlo_error_digital_option_B;

% Display Results
fprintf('\n### MONTE CARLO - BLACK DYNAMICS ###\n')
fprintf('Plain Vanilla Option:\n');
fprintf('Monte Carlo Price (under Black dynamics): %f\n', vanilla_option_B_price);
fprintf('Confidence Interval (under Black dynamics) (%.2f%%): [%.6f, %.6f]\n', confidence_level * 100, lower_bound_plain_vanilla_B, upper_bound_plain_vanilla_B);

fprintf('\nDigital Option:\n');
fprintf('Monte Carlo Price (under Black dynamics): %f\n', digital_option_B_price);
fprintf('Confidence Interval (under Black dynamics) (%.2f%%): [%.6f, %.6f]\n', confidence_level * 100, lower_bound_digital_option_B, upper_bound_digital_option_B);

%% Q.3.4 Overlapping of the confidence intervals ?
% The confidence intervals of level 95% for the price of the plain vanilla 
% option under Local Volatility and Black model overlap and are very
% close. Why do both models give the same result ? Interpretation.
% For the digital option, the confidence intervals don't overlap, they are
% disjoint. Why don't both models give the same result ? Interpretation.
% See slides.













