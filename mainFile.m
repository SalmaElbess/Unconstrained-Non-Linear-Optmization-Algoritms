%%1D minimization 
clc; clear;
f = @(lambda) lambda^5 - 5*(lambda^3) - 20*lambda + 5;
diff_f = @(lambda) 5*lambda^4 - 15*(lambda^2) - 20;
A = 0;
t = 0.5;
eps = 0.0001;
lambda_star_quad = quad_1d_min(f,A,t,eps);
eps1 = 0.05;
eps2 = 0.05;
lambda_star_cub = cub_1d_min(f,diff_f,A,t,eps1,eps2);
[fib_a,fib_b] = fibonacci_1d_minimization(f,6,0,1);
[golden_a,golden_b] = golden_ration_1d_minimization(f,6,0,1);
%% Rosenbrock 
problem_number = 1;
f = @(x) 100*((x(2) - x(1)^2)^2) +(1-x(1))^2;
grad_f = @(x) [-400*(x(2) - x(1)^2)*x(1) - 2*(1-x(1)) ; 200*(x(2)-x(1)^2)];
hessian_f = @(x) [-400*(x(2) - x(1)^2)+ 800*x(1)^2,-400*x(1);-400*x(1),200];
eps = 10^-2; start = [-1.2;1];
%Quasi-Newton Method.
[quasi_i_1,quasi_sol_1,quasi_val_1,quasi_time_1] = quasi_newton(f,grad_f,start,eps,problem_number);
%Fletcher-Reeves CG Method.
[fr_i_1,fr_sol_1,fr_val_1,fr_time_1] = FR(f,grad_f,hessian_f,start);
%Marquardt Method.
[marq_i_1, marq_sol_1,marq_val_1,marq_time_1] = marq(f,grad_f,hessian_f,start,eps,problem_number);
%% powell
problem_number = 2;
f = @(x) (x(1) + 10*x(2))^2 + 5*(x(3)-x(4))^2 +(x(2)-2*x(3))^4 + 10*(x(1)-x(4))^4;
grad_f = @(x) [2*(x(1)+10*x(2)) + 40*(x(1)-x(4))^3 ;
20*(x(1)+10*x(2)) + 4*(x(2)-2*x(3))^3;
10*(x(3)-x(4)) - 8*(x(2)-2*x(3))^3;
-10*(x(3)-x(4)) - 40*(x(1)-x(4))^3];
hessian_f = @(x)[2+120*(x(1)-x(4))^2,20,0,-120*(x(1)-x(4))^2;
    20,200+12*(x(2)-2*x(3))^2,-24*(x(2)-2*x(3))^2,0;
    0,-24*(x(2)-2*x(3))^2,10 + 48*(x(2)-2*x(3))^2,-10;
    -120*(x(1)-x(4))^2, 0, -10, 10+120*(x(1)-x(4))^2];
eps = 10^-4; start = [3;-1;0;1];
%Quasi-Newton Method.
[quasi_i,quasi_sol,quasi_val,quasi_time] = quasi_newton(f,grad_f,start,eps,problem_number);
%Fletcher-Reeves CG Method.
[fr_i,fr_sol,fr_val,fr_time] = FR(f,grad_f,hessian_f,start);
%Marquardt Method.
[marq_i, marq_sol,marq_val,marq_time] = marq(f,grad_f,hessian_f,start,eps,problem_number);