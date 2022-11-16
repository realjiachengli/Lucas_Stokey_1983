function Sol = time0_problem_fminsearch(x0, param)
    % This function solves the time0 social planner's problem where the continuation value
    % function is given by solution in time1_problem.
    % x0 is a struct with all initial values of the model.

    % unpack
    nS = param.nS;
    x_grid = param.x_grid;
    gfunc = param.gfunc;
    opts = param.opts;

    b0 = x0.b0;
    s0 = x0.s0;
    g0 = gfunc(s0);

    % first solve for the continuation bellman equation
    sol = time1_problem(param);

    disp('Time-0 reached.')

    % unpack solution
    Vend = sol.V;
    policy_n = sol.policy_n;
    policy_x_prime = sol.policy_x_prime;

    % define objective2solve s a func of choices for a given state (b0, s0)
    objective2max = @(z) - objective0(z,b0,s0,Vend,param);

    % create bounds for maximization
    nbounds = [0,1];
    xbounds = [min(x_grid), max(x_grid)];
    bounds = [nbounds;repmat(xbounds,[nS-1,1])];
    UB = bounds(:,2);
    LB = bounds(:,1);

    % define initial value as solution of the time1 problem
    n0_init = policy_n(1,s0);
    x_prime0_init = policy_x_prime(:,1,s0);
    z0 = [n0_init;x_prime0_init(1:end-1)];

    % maximize
    [z_sol,fval] = fminsearchbnd(objective2max,z0,LB,UB);

    % get solution as policy
    n0 = z_sol(1);
    x_prime0_ExceptLast = z_sol(2:end);

    [MU_c0,MU_l0,~] = util(n0-g0,1-n0,param);

    x0 = MU_c0*b0;
    x_prime0_Last = find_x_prime_Last(x0,s0,z_sol,MU_c0,MU_l0,param);

    % return everything from time0 and time1 problem
    Sol.V = Vend;
    Sol.W = - fval;
    Sol.policy_n = policy_n;
    Sol.policy_x_prime = policy_x_prime;
    Sol.policy_n0 = n0;
    Sol.policy_x_prime = [x_prime0_ExceptLast;x_prime0_Last];

end


function Val = objective0(z, b0, s0, Vend, param)
    % This function defines the objective of the initial Ramsey planner
    % state variables: b0 and s0;
    % choice variables: z including n and x_prime_ExceptLast;
    % implementability constraint yields x_prime_Last;
    % value function from the continuation ramsey problem: Vend;

    % z is the choice variable as a vector
    n0 = z(1);
    x_prime_ExceptLast = z(2:end);

    beta = param.beta;
    Trans = param.Trans;
    gfunc = param.gfunc;
    x_grid = param.x_grid;
    nS = param.nS;

    g0 = gfunc(s0);   % government spending at initial state s0
    [MU_c0, MU_l0, U] = util(n0-g0,1-n0,param);   % utility


    % get x_prime_Last from implementability constraint
    x0 = MU_c0*b0;
    x_prime_Last = find_x_prime_Last(x0,s0,z,MU_c0,MU_l0,param);
    x_prime = [x_prime_ExceptLast; x_prime_Last];

    % restrict x_prime_Last to the grid range
    if (x_prime_Last<min(x_grid)) || (x_prime_Last>max(x_grid))
        Val = -1e9;
        return;
    end

    % evaluate Vend(x',s')
    Vend_Val = zeros(nS,1);

    for s_prime =1:nS
        x = x_prime(s_prime);
        Vend_Val(s_prime) = interplinear(x,s_prime,Vend,param);
    end

    % finally calculate objective
    Val = U + beta*Trans(s0,:)*Vend_Val;
end



function sol = time1_problem(param)
    % This function solves the continuation Ramsey problem through backward induction
    % and value function iteration until period 1 given parameters set in param.

    nS = param.nS;
    x_grid = param.x_grid;
    MaxIter = param.MaxIter;
    opts = param.opts;
    tol = param.tol;        % stopping threshold

    gfunc = param.gfunc;
    nx = length(x_grid);
    
    Vend = zeros(nx,nS);
    Vcont = zeros(nx,nS);

    % define policy functions
    policy_n = zeros(nx,nS);
    policy_x_prime = zeros(nS,nx,nS);

    % the initial value for objective2solve is n = 0 and x' = mean(x_grid) at steady state for all s';
    % then the initial value for period t iteration is the solution in period t+1; we denote
    % all choice variables z = [n;x_prime_exceptLast]
    z0 = [ ; repmat(mean(x_grid),[nS-1,1])];

    % begin value function iteration
    for iter=1:MaxIter
        for s=1:nS
            for xk=1:length(x_grid)

                % unpack grid point
                x = x_grid(xk);

                % define objective2solve s a func of choices for a given state (x, s)
                objective2min = @(z) - objective(z, x, s, Vend, param);

                % create bounds for maximization
                nbounds = [0,1];
                xbounds = [min(x_grid), max(x_grid)];
                bounds = [nbounds;repmat(xbounds,[nS-1,1])];
                UB = bounds(:,2);
                LB = bounds(:,1);

                % update initial policy as optimum from last iteration
                if iter>1
                    n0 = policy_n(xk,s);
                    x_prime0 = policy_x_prime(:,xk,s);
                    z0 = [n0;x_prime0(1:end-1)];
                end

                % maximize
                [z_sol,fval] = fminsearchbnd(objective2min,z0,LB,UB);

                Vcont(xk,s) = -fval;

                n = z_sol(1);
                g = gfunc(s);
                [MU_c, MU_l, ~] = util(n-g,1-n,param);   % utility

                x_prime_Last = find_x_prime_Last(x,s,z_sol,MU_c,MU_l,param);

                % pack policy function
                policy_n(xk,s) = z_sol(1);
                policy_x_prime(:,xk,s) = [z_sol(2:end);x_prime_Last];

            end
        end
        % report progress
        dist = max(max(abs(Vcont-Vend)));
        fprintf('Time-1 progress: iter=%d, dist=%f\n', iter, dist)

        % stopping criterion
        if dist < tol
            break;
        end

        Vend = Vcont;
    end

    sol.V = Vcont;
    sol.policy_n = policy_n;
    sol.policy_x_prime = policy_x_prime;
end


function Val = objective(z, x, s, Vend, param)
    % This function defines the objective of the continuation Ramsey planner
    % state variables: x and s;
    % choice variables: z including n and x_prime_ExceptLast;
    % implementability constraint yields x_prime_Last;
    % value function from the last iteration: Vend;

    % z is the choice variable as a vector
    n = z(1);
    x_prime_ExceptLast = z(2:end);


    beta = param.beta;
    Trans = param.Trans;
    gfunc = param.gfunc;
    x_grid = param.x_grid;
    nS = param.nS;

    g = gfunc(s);   % government spending at state s
    [MU_c, MU_l, U] = util(n-g, 1-n,param);   % utility

    % get x_prime_Last from implementability constraint
    x_prime_Last = find_x_prime_Last(x,s,z,MU_c,MU_l,param);
    x_prime = [x_prime_ExceptLast; x_prime_Last];

    % restrict x_prime_Last to the grid range
    if x_prime_Last<min(x_grid) | x_prime_Last>max(x_grid)
        Val = -1e9;
        return;
    end

    % evaluate Vend(x',s')
    Vend_Val = zeros(nS,1);
    for s_prime =1:nS
        x = x_prime(s_prime);
        Vend_Val(s_prime) = interplinear(x,s_prime,Vend,param);
    end

    % finally calculate objective
    Val = U + beta*Trans(s,:)*Vend_Val;
end


function x_prime_Last = find_x_prime_Last(x,s,z,MU_c,MU_l,param)
    % calculates the last element of x' from the implementability constraint
    Trans = param.Trans;
    beta = param.beta;

    % unpack z
    n = z(1);
    x_prime_ExceptLast = z(2:end);

    % subtract RHS from LHS
    LminusR = x - MU_c + MU_l*n - beta*Trans(s,1:end-1)*x_prime_ExceptLast;
    x_prime_Last = LminusR/(beta*Trans(s,end));
end


function [MU_c, MU_l, U] = util(c, l, param)
    % compute utility, marginal utility of consumption and marginal utility of labor
    gamma = param.gamma;
    if c<0
        U = -1e9;
    else
        U = log(c) + gamma*log(l);
    end
    MU_c = 1/c;
    MU_l = gamma/l;
end


function Val = interplinear(x,s_prime,Vend,param)
    % perform a linear interpolation
    Vend_s = Vend(:,s_prime);

    x_grid = param.x_grid;
    i = lookup(x_grid,x,3);
    weightL = (x-x_grid(i)) / (x_grid(i+1)-x_grid(i));
    % evaluate Val as weighted sum of two endpoints
    Val = (1-weightL)*Vend_s(i) + weightL*Vend_s(i+1);
end