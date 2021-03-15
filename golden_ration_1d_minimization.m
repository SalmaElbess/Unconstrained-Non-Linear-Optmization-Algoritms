function[a,b] = golden_ration_1d_minimization(f,n,a,b)
N = n+1;
for i = N:-1 :3
L2 = 0.382*(b-a);
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

elseif f1 <f2
        b = x2;

elseif f1==f2
    a = x1;
    b=x2;
end
end
end