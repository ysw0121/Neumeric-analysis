syms x f;
f(x)=2*log(x)+log(3);
g(x)=x-(f(x)-x)^2/(f(f(x))-2*f(x)+x);
eps=10^(-10);
k=0;
step=20;
a=3;b=4;
x0=a;
x1=g(x0);
if(x1>4)
    fprintf('在[3,4]不收敛\n');
   return;
end
fprintf('%d   %.8f    %.8f\n',k,x0,x1);
while(abs(x1-x0)>eps&&step>0)
    x1=g(x0);
    x0=x1;
    k=k+1;
    step=step-1;
    fprintf('%d   %.8f    %.8f\n',k,x0,x1);
end
if(abs(x1-x0)<=eps)
    fprintf('结果是%.8f\n',x1);
end
if(step==0)
    fprintf('步数不足\n');
end