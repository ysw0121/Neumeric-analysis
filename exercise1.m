%slCharacterEncoding('UTF-8')
% clear;clc
% disp('输入矩阵：')
% A=input('矩阵A=')
% B=A';
% C=inv(A);
% B,C


%   二分法求函数的实根

% % clear;
% % clc;
% format short
% s=input('请输入函数表达式：f = ','s');
% %f=inline(s);
% f=eval(['@(x)',s]);
% a = input('请输入区间左端点:a=');
% b = input('请输入区间右端点:b=');
% if a>b
%     fprintf('区间表示不正确，请重新开始');
%     return;
% end
% if f(a)*f(b)>0
%     fprintf('两端点不存在变号零点，重新开始');
%     return;
% end
% eps=input('请输入停止精度要求:eps=');  %%“|b-x|<=eps”
% k=1;
% x=(a+b)/2;
% fprintf('        k        a        f(a)       b        f(b)       xk        f(xk)\n ');
% T=[k,a,f(a),b,f(b),x,f(x)];
% %while abs(T(k,4)-T(k,6))>eps/2
% while abs(T(k,2)-T(k,4))>eps
%     k=k+1;
%     if  f(x)*f(a)==0
%         a=a;
%         b=x;
%         x=(a+b)/2;
%         T=[T;k,a,f(a),b,f(b),x,f(x)];
%         break
%     elseif  f(x)*f(a)>0
%         a=x;
%         b=b;
%         x=(a+b)/2;
%         T=[T;k,a,f(a),b,f(b),x,f(x)];
%     elseif  f(x)*f(a)<0
%         a=a;
%         b=x;
%         x=(a+b)/2;
%         T=[T;k,a,f(a),b,f(b),x,f(x)];
%     end
% 
% end
% disp(T);
% fprintf('经过%d次迭代，函数方程根的近似解为：x=%.8f\n',k-1,T(k-1,6))


% %迭代法求解x=f(x)

% clc;
% clear;
% syms x;
% f=input("请输入迭代方程(自变量为x,如1/3*(x^3+1)):  ");
% p0=input("请输入不动点迭代法的初始值：  ");
% perror=input("请输入允许的误差值：  ");
% maxK=input("请输入最大迭代次数：  ");
% 
% [p,k,Y]=FPM(f,p0,perror,maxK);
% 
% DP=sprintf("使用不动点迭代法迭代%d次，计算%s=x在%g附近的解为：%g",k,f,p0,p);
% disp(DP);
% fprintf("迭代值如下：");
% disp(Y);
% 
% function [p,k,Y]=FPM(f,p0,perror,maxK)
% %p0表示迭代初始值
% %f表示迭代公式函数
% %maxK表示规定的最大迭代次数
% %pererr表示允许误差
% %k表示最终迭代的次数
% %p表示最终迭代的值
% %Y用来记录每次迭代过程的迭代值
%     syms x;
%     P(1)=p0;
%     k=2;
%     P(k)=subs(f,x,P(k-1));      %迭代
%     while k<=maxK
%         err=abs(P(k)-P(k-1));    %err表示相邻的迭代值的差值
%         if(err<perror)
%             fprintf('迭代%d次即可满足允许误差值退出\n',k-1);
%             break;
%         end
%         k=k+1;
%         P(k)=subs(f,x,P(k-1));
%     end         %共迭代了k-1次
%     if(k-1==maxK) 
%         disp("超过最大迭代次数！");
%     end
%     p=P(k); 
%     k=k-1;
%     Y=P;
% end

%Newton 迭代法求根



