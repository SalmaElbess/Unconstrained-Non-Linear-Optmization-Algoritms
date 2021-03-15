function [i, x_old,f_old,imp_time] = quasi_newton(f,grad_f,x_new,eps,problem_number)
%intial values
B = eye(length(x_new));
i = 0;
imp_start = now;
%start algorithm
while (true)
    i = i+1;
    
    %update old values
    x_old = x_new;
    f_old = f(x_old);
    grad_old = grad_f(x_old);
    
    %stoping criteria
    if norm(grad_old) < eps
        break;
    end
    
    %get S
    S =  - B\grad_old;
    
    %get lambda
    if problem_number ==1
        f_lambda = @(lambda) 100*(((x_old(2)+lambda*S(2)) - (x_old(1)+lambda*S(1))^2)^2) +(1-(x_old(1)+lambda*S(1)))^2;
    else
        f_lambda = @(lambda) ((x_old(1)+lambda*S(1)) + 10*(x_old(2)+lambda*S(2)))^2 + 5*((x_old(3)+lambda*S(3))-(x_old(4)+lambda*S(4)))^2 +((x_old(2)+lambda*S(2))-2*(x_old(3)+lambda*S(3)))^4 + 10*((x_old(1)+lambda*S(1))-(x_old(4)+lambda*S(4)))^4;
    end
    [lower_limit,upper_limit] = fibonacci_1d_minimization(f_lambda,11,0,1);
    lambda_star = (upper_limit+lower_limit)/2;
    
    %get new values
    x_new = x_old + lambda_star.*(S);
    grad_new = grad_f(x_new);
    
    %update B
    y = grad_new - grad_old;
    g=x_new - x_old;
    B = B + (y-B*g)*(y-B*g)'/((y-B*g)'*g);
end
%end of algoritm
imp_end= now;
imp_time = (imp_end - imp_start)*86400;
end
