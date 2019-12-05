function  [x_optimal cash_optimal] = strat_min_variance(x_init, cash_init, mu, Q, cur_prices)

    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
  
    n = length(x_init);
    x_opt = zeros(n,1);
    
    % Optimization problem data
    lb = zeros(n,1);
    ub = inf*ones(n,1);
    A  = ones(1,n);
    b  = 1;
    
    % Compute minimum variance portfolio
    cp_v = Cplex('min_var');
    cp_v.addCols(zeros(n,1), [], lb, ub);
    cp_v.addRows(b, A, b);
    cp_v.Model.Q = 2*Q;
    cp_v.Param.qpmethod.Cur = 6; % concurrent algorithm
    cp_v.Param.barrier.crossover.Cur = 1; % enable crossover
    cp_v.DisplayFunc = []; % disable output to screen
    cp_v.solve();
    
    %optimized weight from solver
    w_opt = cp_v.Solution.x;
    
    %profolio value
    pro_v = dot(x_init,cur_prices);

    for i = 1:n
        x_opt(i) = floor(w_opt(i) * pro_v / cur_prices(i)); 
    end
    
    x_opt(x_opt<0)=0;
    
    %calculate total stock sold/bought
    x_delta = x_init-x_opt;
    %calcualte transaction fee
    tran_fee = dot(abs(x_delta),cur_prices)*0.005;
    %rebalance the cash
    cash_opt = cash_init + dot(x_delta, cur_prices) - tran_fee;
    
    %buy more stocks with leftover cash
    while cash_opt>0
        price_delta = abs(cur_prices-abs(cash_opt));
        [max_diff,index]=max(price_delta);
        x_opt(index) = x_opt(index) +1;
        cash_opt = cash_opt - dot(cur_prices(index),1.005);
    end
    
    %check for negative balance and sold stocks bought accordingly
    while cash_opt <0
        
         %calculate the diff btw cur_price and cash
        price_delta = abs(cur_prices-abs(cash_opt));
        
        %1 == stocks bought, inf == stocks sold
        x_delta = x_init-x_opt;
        x_check=zeros(size(x_delta));
        for i=1:length(x_delta)
            if x_delta(i)<0
                x_check(i)=1;
            else
                x_check(i)=inf;
            end     
        end
        %find the closest stock price
        x_diff = x_check'.*price_delta;
        [min_diff,index]=min(x_diff);
        
        %rebalance the cash and holdings
        x_opt(index) = x_opt(index) -1;
        cash_opt = cash_opt + dot(cur_prices(index),1.005);

    end
    
    %buy more stocks with leftover cash
    while cash_opt>0
        price_delta = abs(cur_prices-abs(cash_opt));
        [max_diff,index]=max(price_delta);
        cash_tmp = cash_opt - dot(cur_prices(index),1.005);
        
        if cash_tmp >=0
            x_opt(index) = x_opt(index) +1;
            cash_opt = cash_tmp;
        else
            break;
        end
    end
    
    %disp(cash_opt);
    x_optimal = x_opt;
    cash_optimal = cash_opt;

   
    
end

