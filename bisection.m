
%clc;clear all;
s=input('请输入函数表达式：f = ','s');
f=inline(s);
a = input('请输入区间左端点:a=');
b = input('请输入区间右端点:b=');
eps=input('请输入停止精度要求:eps='); 
k=1;
x=(a+b)/2;
fprintf('k                    a                     f(a)                    b                f(b)                      xk                        f(xk)\n ');
if(f(a)*f(b)>0)
    disp('不符要求，失败');
    return;
end
x=(a+b)/2;
while(abs(a-b)>=eps)   
    fprintf('%d       %.8f       %.8f          %.8f           %.8f           %.8f           %.8f\n',k,a,f(a),b,f(b),x,f(x));
    x=(a+b)/2;
    if(f(a)*f(x)<0)
        a=a;
        b=x;
    elseif(f(b)*f(x)<0)
            b=b;
            a=x;
    end
    k=k+1;
end
fprintf('结果是%.8f\n',x);