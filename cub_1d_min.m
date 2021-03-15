function lambda_star = cub_1d_min(f,diff_f,A,t,eps1,eps2)
diff_t = diff_f(t);
while diff_t < 0
    t = t*2;
    diff_t = diff_f(t);
end
B = t;
i = 0; % number of iterations
while true
    i = i+1;
    fa = f(A); fb = f(B); diff_fa = diff_f(A); diff_fb = diff_f(B);
    Z = (3*(fa- fb)/(B-A)) +diff_fa +diff_fb;
    Q = sqrt(Z^2 - (diff_fa*diff_fb));
    lambda_star1 = A + ((diff_fa+Z+Q)*(B-A)/(diff_fa+diff_fb + 2*Z));
    lambda_star2 = A + ((diff_fa+Z-Q)*(B-A)/(diff_fa+diff_fb + 2*Z));
    if lambda_star1 <= B && lambda_star1 >= A
        lambda_star = lambda_star1;
    else
        lambda_star = lambda_star2;
    end
    b = (diff_fa*(B^2)+diff_fb*(A^2) + 2*A*B*Z)/((A-B)^2);
    c = -1*(Z*(A+B) + B*diff_fa + A*diff_fb)/((A-B)^2);
    d = (2*Z + diff_fa +diff_fb)/(3*(A-B)^2);
    a = fa - b*A - c*(A^2) - d*(A^3);
    h = a + b*lambda_star + c*lambda_star^2 + d*lambda_star^3;
    f_lambda = f(lambda_star);
    diff_f_lambda = diff_f(lambda_star);
    if abs(diff_f_lambda) <= eps2 && abs((h-f_lambda)/f_lambda) <= eps1
        break;
    end
    if diff_f_lambda > 0
        B = lambda_star;
    else
        A = lambda_star;
    end
end
end