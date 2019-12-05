function  [x_optimal cash_optimal] = strat_robust_optim(x_init, cash_init, mu, Q, cur_prices)

    % Add PATH to CPLEX and ipopT solvers
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    
    disp(cur_prices)

    n = length(x_init);
    
    % Define initial portfolio ("equally weighted" or "1/n portfolio")
    w0 = ones(n,1) ./ n;

    % Bounds on variables
    
    lb_rMV = zeros(n,1);
    ub_rMV = inf*ones(n,1);

    % Required portfolio robustness
    var_matr = diag(diag(Q));
    % Target portfolio return estimation error is return estimation error of 1/n portfolio
    rob_init = w0' * var_matr * w0; % return estimation error of initial portfolio
    rob_bnd = 0.01*rob_init; % target return estimation error..make this small

    % Compute minimum variance portfolio
    cplex_minVar = Cplex('MinVar');
    
    cplex_minVar.DisplayFunc = [];
    
    cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
    cplex_minVar.addRows(1, ones(1,n), 1);
    cplex_minVar.Model.Q = 2*Q;
    cplex_minVar.Param.qpmethod.Cur = 6;
    cplex_minVar.solve();
    cplex_minVar.Solution;
    
    w_minVar = cplex_minVar.Solution.x; % asset weights
    
    ret_minVar = dot(mu, w_minVar);
    
    % Target portfolio return is return of minimum variance portfolio
    Portf_Retn = 10*ret_minVar; %make this greater for profit maxmization
    
    % Objective function
    f_rMV  = zeros(n,1);
    % Constraints
    A_rMV  = sparse([  mu';...
                 ones(1,n)]);
    lhs_rMV = [Portf_Retn; 1];
    rhs_rMV = [inf; 1];
    % Initialize CPLEX environment
    cplex_rMV = Cplex('Robust_MV');
    
    cplex_rMV.DisplayFunc = [];
    
    % Add objective function and variable bounds
    cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
    % Add constraints
    cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
    % Add quadratic objective
    cplex_rMV.Model.Q = 2*Q;
    % Add quadratic constraint on return estimation error (robustness constraint)
    Qq_rMV = var_matr;
    cplex_rMV.addQCs(zeros(size(f_rMV)), Qq_rMV, 'L', rob_bnd, {'qc_robust'});
    % Set CPLEX parameters
    cplex_rMV.Param.threads.Cur = 4;
    cplex_rMV.Param.timelimit.Cur = 60;
    cplex_rMV.Param.barrier.qcpconvergetol.Cur = 1e-12; % solution tolerance
    
    
    
    cplex_rMV.solve();   
    cplex_rMV.Solution;
    
    
    w_rMV = cplex_rMV.Solution.x;

    
    %profolio value = total holding - interest of loan
    pro_v = dot(x_init,cur_prices);
    
    x_opt = zeros(n,1);
    
    for i = 1:n
        x_opt(i) = floor(w_rMV(i) * pro_v / cur_prices(i)); 
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

