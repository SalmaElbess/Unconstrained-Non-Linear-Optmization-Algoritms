function [i, x_old,f_old,imp_time] = marq(f,grad_f,hessian_f,x_new,eps,problem_number)
%initial values
imp_start = now;
alpha1 = 10^4; c1 = 0.5; c2 = 2;
i = 0;
%start algorithm
while (true)
    i = i+1;
    %update old values
    x_old = x_new; f_old = f(x_old);
    grad = grad_f(x_old); f_hess = hessian_f(x_old);
    %stoping criteria
    if norm(grad) < eps
        break;
    end
    
    %get S
    I = eye(size(f_hess));
    S = -1*(f_hess+ alpha1*I)^-1 * grad;
    
    %get lambda
    if problem_number ==1
        f_lambda = @(lambda) 100*(((x_old(2)+lambda*S(2)) - (x_old(1)+lambda*S(1))^2)^2) +(1-(x_old(1)+lambda*S(1)))^2;
    else
        f_lambda = @(lambda) ((x_old(1)+lambda*S(1)) + 10*(x_old(2)+lambda*S(2)))^2 + 5*((x_old(3)+lambda*S(3))-(x_old(4)+lambda*S(4)))^2 +((x_old(2)+lambda*S(2))-2*(x_old(3)+lambda*S(3)))^4 + 10*((x_old(1)+lambda*S(1))-(x_old(4)+lambda*S(4)))^4;
    end
    
    [lower_limit,upper_limit] = fibonacci_1d_minimization(f_lambda,11,0,1);
    lambda_star = (upper_limit+lower_limit)/2;
    
    %update new values
    x_new = x_old + lambda_star*S;
    f_new = f(x_new);
    if f_new < f_old
        alpha1 = alpha1*c1;
    else
        alpha1 = c2*alpha1;
    end
end
imp_end= now;
imp_time = (imp_end - imp_start)*86400;
end