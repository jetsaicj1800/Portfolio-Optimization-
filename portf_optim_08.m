clc;
clear all;
format long

global Q

% Input files
%input_file_prices  = 'Daily_closing_prices.csv';
input_file_prices  = 'Daily_closing_prices20082009.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2007)); %2014 or 2007
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;
% Annual risk-free rate for years 2008-2009 is 4.5%
r_rf2008_2009 = 0.045;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe' 'strat_equal_risk_contr' 'strat_lever_equal_risk_contr' 'strat_robust_optim'};
%strategy_functions = {'strat_equal_risk_contr','strat_lever_equal_risk_contr'};

strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Minimum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio' 'Equal Risk Contributions Portfolio' 'Leveraged Equal Risk Contributions Portfolio' 'Robust Optimization Portfolio'};
%N_strat = 6; % comment this in your code
N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);

for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==08)
       cur_year  = 08 + floor(period/7);
   else
       cur_year  = 2008 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:N_strat)

      % Get current portfolio positions
      if(period==1)
         %double the initial holding for leverage equal risk
         if strategy == 6
             curr_positions = 2*curr_positions;
         else
             curr_positions = init_positions;
         end
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end

      % Compute strategy
      [x{strategy,period} cash{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices);

      % Verify that strategy is feasible (you have enough budget to re-balance portfolio)
      % Check that cash account is >= 0
      % Check that we can buy new portfolio subject to transaction costs

      %%%%%%%%%%% Insert your code here %%%%%%%%%%%%
      
      if cash{strategy,period} < 0
          error(' we can not have negative balance in cash');
      end

      % Compute portfolio value
      portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};
      
      %subtract the loan from the portfolio value calculation
      if strategy == 6
          portf_value{strategy}(day_ind_start:day_ind_end) = portf_value{strategy}(day_ind_start:day_ind_end) - init_value;
      end
      
      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end));

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);
   
end

% Plot results
% figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
%plot strategy charts
figure(1);
plot(portf_value{1}(:,:),'r');
hold on;
plot(portf_value{2}(:,:),'g');
hold on;
plot(portf_value{3}(:,:),'c');
hold on;
plot(portf_value{4}(:,:),'b');
hold on;
plot(portf_value{5}(:,:),'m');
hold on;
plot(portf_value{6}(:,:),'k');
hold on;
plot(portf_value{7}(:,:),'y');
hold on;

xlabel('Days');
ylabel('Profolio Value');
title('Performance of Strategies');
legend(strategy_names,'Location','southoutside');

%plot dynamic changes

w_rob = zeros(20,N_periods);
w_var = zeros(20,N_periods);
w_sharp = zeros(20,N_periods);
%change for min-variance
for period = 1:N_periods
    
   if(dates_array(1,1)==08)
       cur_year  = 08 + floor(period/7);
   else
       cur_year  = 2008 + floor(period/7);
   end
   
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   
   cur_prices = data_prices(day_ind_start,:);
   
   w_rob(:,period) = (cur_prices' .* x{7,period})/(cur_prices *x{7,period}) ;
   w_var(:,period) = (cur_prices' .* x{3,period})/(cur_prices *x{3,period}) ;
   w_sharp(:,period) = (cur_prices' .* x{4,period})/(cur_prices *x{4,period});
   
end

figure(2)
plot(w_rob')
xlabel('Period');
ylabel('Weight');
title('Robust Optimization Portfolio Dynamic Change');
legend('Location','eastoutside');

figure(3)
plot(w_var')

xlabel('Period');
ylabel('Weight');
title('Min Variance Strategy Dynamic Change');
legend('Location','eastoutside');

figure(4)
plot(w_sharp')
xlabel('Period');
ylabel('Weight');
title('Max Sharpe Ratio Dynamic Change');
legend('Location','eastoutside');
