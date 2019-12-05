function  [x_optimal cash_optimal] = strat_max_Sharpe(x_init, cash_init, mu, Q, cur_prices)

    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    n = length(x_init);
    x_opt = zeros(n,1);
    
    r_rf = 0.025;
    %r_rf = 0.045;
    %r_rf2008_2009 = 0.045;
    
    newQ = zeros(size(Q)+1);
    newQ(1:size(Q,1),1:size(Q,2)) = Q;
    
    
    % Optimization problem data
    lb = zeros(n+1,1);
    ub = inf*ones(n+1,1);
    A = [[(  mu'- r_rf/252) 0]; [ones(1,n) -1]];
    b = [1; 0]; 
    
    
    %check if A contains negative values
    if all(A(1,:) <=0)
        x_optimal = x_init;
        cash_optimal = cash_init;
        return 
    end
    

   
    % Compute minimum variance portfolio
    cp_s = Cplex('max_sharpe');
    cp_s.Model.sense = 'minimize';
    cp_s.addCols(zeros(n+1,1), [], lb, ub);
    cp_s.addRows(b, A, b);
    cp_s.Model.Q = 2*newQ;
    cp_s.Param.qpmethod.Cur = 6; % concurrent algorithm
    cp_s.Param.barrier.crossover.Cur = 1; % enable crossover
    cp_s.DisplayFunc = []; % disable output to screen
    cp_s.solve();


    y = cp_s.Solution.x(1:n);
    k = cp_s.Solution.x(n+1);
    w_s = y/k;
    
    pro_v = dot(x_init,cur_prices);

    for i = 1:n
        x_opt(i) = floor(w_s(i) * pro_v / cur_prices(i)); 
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