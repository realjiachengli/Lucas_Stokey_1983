function xSS = find_xSS(param)
    % This function find the deterministic steady state of the time1 bellman equation

    beta = param.beta;
    gamma = param.gamma;
    gfunc = param.gfunc;

    % set g to the average of all states
    g = mean(gfunc);

    



end

% function objective


function x_prime_Last = find_x_prime_Last(x,s,z,param)
    % calculates the last element of x' from the implementability constraint
    Trans = param.Trans;
    beta = param.beta;
    g = gfunc(s);   % government spending at state s

    % unpack z
    n = z(1);
    x_prime_ExceptLast = z(2:end);

    [MU_c, MU_l, ~] = util(n-g, 1-n);   % utility

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
