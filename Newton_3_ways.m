syms x f;
f(x)=x^4-4*x^2+4;
x0=1.5;
eps=1e-6;
step=5;
df=diff(f,1);
d=double(subs(df,x0));

fprintf('牛顿法如下：\n');
for i=1:step
    x1=x0-(f(x0)./d);
    df=diff(f,1);
    d=double(subs(df,x0));
    fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x0, f(x0));
     if abs(x1-x0)<eps, break, end
     x0=x1;
end
fprintf('\n');

x0=1.5;
df=diff(f,1);
d=double(subs(df,x0));
fprintf('已知二重根的牛顿法如下：\n');
for i=1:step
    x1=x0-2*(f(x0)./d);
    df=diff(f,1);
    d=double(subs(df,x0));
    fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x0, f(x0));
     if abs(x1-x0)<eps, break, end
     x0=x1;
end
fprintf('\n');


x0=1.5;
df=diff(f,1);
d=double(subs(df,x0));
fprintf('不知二重根的牛顿法如下：\n');
for i=1:step
    d2f=diff(f,2);
    d2=double(subs(d2f,x0));
    df=diff(f,1);
    d=double(subs(df,x0));
    x1=x0-(f(x0)*d2)/(d^2-f(x0)*d2);
    fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x0, f(x0));
    if abs(x1-x0)<eps, break, end
     x0=x1;
end
fprintf('\n');
