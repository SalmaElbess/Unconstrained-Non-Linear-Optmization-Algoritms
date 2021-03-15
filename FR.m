function [i, x_new,f_new,imp_time] = FR(f,grad_f,hessian_f,x_old)
%initial values
imp_start = now;
grad_old = grad_f(x_old);
S = - grad_old/norm(grad_old);
i = 0;

%start algorithm
while true
    i = i+1;
    
    %get lambda
    A = hessian_f(x_old);
    lambda_star = norm(grad_old/norm(grad_old))^2 /(S'*A*S);
    
    %update x
    x_new = x_old + lambda_star.*S;
    grad_new = grad_f(x_new);
    
    %stoping criteria
    if norm(grad_new) < 1e-3
        break;
    end
    
    %update s
    beta = ((norm(grad_new))^2)/ ((norm(grad_old))^2);
    S = -1*grad_new + beta*S;
    S = S/norm(S);
    
    %change inital values
    x_old = x_new;
    grad_old = grad_new;
    
end
%end of algoritm

%get final values
f_new = f(x_new);
imp_end= now;
imp_time = (imp_end - imp_start)*86400;
end