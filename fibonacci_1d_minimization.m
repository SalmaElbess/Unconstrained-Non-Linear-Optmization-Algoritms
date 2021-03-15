function[a,b] = fibonacci_1d_minimization(f,n,a,b)
N = n+1;
%get fibonacci sequence
fib_seq = zeros(1,n+1);
fib_seq(1) = 1; fib_seq(2) = 1;
for i = 3:length(fib_seq)
    fib_seq(i) = fib_seq(i-1) + fib_seq(i-2);
end

for i = N:-1 :3
L2 = fib_seq(i-2)*(b-a)/fib_seq(i);
L1 = b-a;
if L2 > L1/2
    x1 = b-L2;
    x2 = a+L2;
else
    x2 = b-L2;
    x1 = a+L2;
end
f1 = double(f(x1));
f2 = double(f(x2));
if f1 > f2
        a = x1;

else
        b = x2;


end
end
end