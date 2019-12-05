function  [x_optimal cash_optimal] = strat_equal_risk_contr(x_init, cash_init, mu, Q, cur_prices)

    global A_ineq A_eq
    % Add PATH to CPLEX and ipopT solvers
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    addpath('C:\Users\jetsa\Desktop\MIE1622\IPOPT');
    addpath('C:\Users\jetsa\Desktop\MIE1622\Assign2\equal_risk_contrib');
    
    n = length(x_init);
    x_opt = zeros(n,1);
    
    % Equality constraints
    A_eq = ones(1,n);
    b_eq = 1;

     % Inequality constraints
    A_ineq = [];
    b_ineql = [];
    b_inequ = [];
           
    % Define initial portfolio ("equally weighted" or "1/n portfolio")
    w0 = repmat(1.0/n, n, 1);

    options.lb = zeros(1,n);       % lower bounds on variables
    options.lu = ones (1,n);       % upper bounds on variables
    options.cl = [b_eq' b_ineql']; % lower bounds on constraints
    options.cu = [b_eq' b_inequ']; % upper bounds on constraints

    % Set the IPOPT options
    options.ipopt.jac_c_constant        = 'yes';
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.tol                   = 1e-10;
    options.ipopt.print_level = 0;
    
    %options.ipopt.derivative_test       = 'first-order';
    %options.ipopt.derivative_test_print_all = 'yes';

    % The callback functions
    funcs.objective         = @computeObjERC;
    funcs.constraints       = @computeConstraints;
    funcs.gradient          = @computeGradERC;
    funcs.jacobian          = @computeJacobian;
    funcs.jacobianstructure = @computeJacobian;
    
    % Run IPOPT
    [w_erc info] = ipopt(w0',funcs,options);

    %profolio value
    pro_v = dot(x_init,cur_prices);

    for i = 1:n
        x_opt(i) = floor(w_erc(i) * pro_v / cur_prices(i)); 
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