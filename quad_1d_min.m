function lambda_star = quad_1d_min(f,A,t,eps)
fa= f(A);
f1 = f(t);
if f1 > fa
    fc = f1;
    fb = f(t/2);
    t = t/2;
elseif f1<fa
    while(true)
        fb = f1;
        f2 = f(2*t);
        if f2 > f1
            fc = f2;
            break;
        else
            f1 = f2;
            t = 2*t;
        end
    end
end
B = t;
C = 2*t;
if A == 0
    lambda_star = (4*fb - 3*fa - fc)*t/(4*fb - 2*fc - 2*fa);
else
    
    lambda_star = fb*(C^2 - A^2) + fa*(B^2 - C^2) + fc*(A^2 - B^2)/2*(fa*(B-C) +fb*(C-A) + fc*(A-B));
end
syms a b c
[sola,solb, solc] = solve(fa == a+b*A + c*(A^2),fb == a+b*B + c*(B^2),fc == a+b*C + c*(C^2),[a b c]);
h = double(sola)+double(solb)*lambda_star + double(solc) *(lambda_star^2);
f_lambda = f(lambda_star);
i = 1;
%%refitting
if abs((h - f_lambda)/f_lambda) > eps
    while abs((h - f_lambda)/f_lambda) > eps
        i = i+1;
    if lambda_star > B
        if f_lambda < fb %neglect old A
            A = B;
            B = lambda_star;
        else %neglect old C
            C = lambda_star;
        end
    else
       if f_lambda < fb %neglect old C
           C = B;
            B = lambda_star;
        else %neglect old A
            A = lambda_star;
        end
    end
    fb = f(B); fc = f(C); fa = f(A);
    syms a b c
[sola,solb, solc] = solve(fa == a+b*A + c*(A^2),fb == a+b*B + c*(B^2),fc == a+b*C + c*(C^2),[a b c]);
lambda_star = -1*double(solb)/(2*double(solc));
    h = double(sola)+double(solb)*lambda_star + double(solc) *(lambda_star^2);
    f_lambda = f(lambda_star);
    n = abs((h - f_lambda)/f_lambda);
    end
end
end