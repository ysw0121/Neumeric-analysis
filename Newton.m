syms x f;
f(x)=(2-x)*exp(2*x-1)+1;
x0=2;
eps=10^(-22);
k=0;
step=20;

df=diff(f,1);
pretty(df);
d=double(subs(df,x0));

disp('k         xk        |x1-x0|    ');

x1=x0-(f(x0)./d);
if(abs(x1-x0)<eps)
    fprintf('只有一次');
    fprintf('%d  %f   %f\n',k,x1,abs(x1-x0));
    return
end

fprintf('%d  %f   %f\n',k,x1,abs(x1-x0));
while(abs(x1-x0)>=eps&&step>0)
    x1=x0-(f(x0)./d);
    delta=abs(x1-x0);
    x0=x1;
    df=diff(f);
    d=double(subs(df,x0));
    k=k+1;
    step=step-1;
    fprintf('%d  %.8f   %.8f\n',k,x1,delta);
end
if(abs(x1-x0)<eps)
    %fprintf('%d  %.8f   %.8f\n',k,x1,delta);
    fprintf('近似解是%.8f\n',x1);

elseif(step==0)
    fprintf('此精度和步数下不可求出值\n');
end
