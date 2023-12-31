
%切比雪夫多项式
% syms t x;
% n=6;
% c={'r','b','k','g','c','y','m'};
% subplot(1,2,1);
% for i=0:n
%     %t=chebyshevT(i,x);
%     t = cos(i*acos(x));
%     h=fplot(t,[-1 1],'LineWidth',1);
%     set(h,'Color',c{i+1});
%     hold on;
% end
% legend('n=0','n=1','n=2','n=3','n=4','n=5','n=6','Location','NorthEastOutside')
% title('Chebyshev polynomials');
% xlabel('x');

% 勒让德多项式
% subplot(1,2,2);
% for i=0:n
%     %t=legendreP(i,x);
%     t= 1/(2^i*prod(1:i))*diff((x^2-1)^i,x,i);
%     h=fplot(t,[-1 1],'LineWidth',1);
%     set(h,'Color',c{i+1});
%     hold on;
% end
% legend('n=0','n=1','n=2','n=3','n=4','n=5','n=6','Location','NorthEastOutside')
% title('Legendre polynomials');
% xlabel('x');

%%---------------------------------------------------------------------------------------------------------------

%拉格朗日插值
% function f=lagrange(x,x0,y0)
% f=0;
% n=length(x0);
% for i=1:n
%     t=1;
%     for j=1:n
%         if j~=i
%             t=t*(x-x0(j))/(x0(i)-x0(j));
%         end
%     end
% f=f+y0(i)*t;
% end
% end

%%---------------------------------------------------------------------------------------------------------------

%二分法求解方程
% syms x;
% k=1;
% a=1;b=2;
% y=@(x)x^3-x-1;
% eps=1e-5;
% x=(a+b)/2;
% while(abs(a-b)>=eps)
%      fprintf('%d       %.8f       %.8f          %.8f           %.8f           %.8f           %.8f\n',k,a,y(a),b,y(b),x,y(x));
%      x=(a+b)/2;
%      if(y(x)*y(a)<0)
%          b=x;
%     elseif(y(x)*y(b)<0)
%         a=x;
%      end
%      k=k+1;
% end
% fprintf('result is %.8f',x);
% 
% % function y0=y(x)
% % y0=x^3-x-1;
% % end

%%---------------------------------------------------------------------------------------------------------------

%不动点解方程
% syms x;
% f(x)=@(x)x^3-x-1;
% g(x)=(x+1)^(1/3);
% step=10;
% fvalue=zeros(step,1);
% x=1.5;
% for i=1:step
%     x=g(x);
%     fvalue(i)=abs(f(x));
%     fprintf("x=%.4f f=%.4f\n",x,abs(f(x)));
% end
% fvalue=log(fvalue);
% plot(1:step, fvalue, 'ob-')
% xlabel('k'); ylabel('|f(x)|');
% title('函数值下降曲线');

%%---------------------------------------------------------------------------------------------------------------

% Aitken加速迭代
% clear all;
% close all;
% 
% f = @(x) x^3-x-1;      % 原函数
% g = @(x) (x+1)^(1/3);  % 迭代函数 
% maxit = 10;            % 最大迭代步数
% 
% fprintf('迭代函数: g=%s\n', char(g));
% 
% fvalue = zeros(maxit,1);
% 
% x = 1.5; % 取中点为迭代初始值
% for k = 1 : maxit
%     x = g(x);
%     fvalue(k) = abs(f(x));
%     fprintf('k=%d, x=%.4f, f(x)=%.2e\n',k, x, fvalue(k));
% end
% 
% % 绘图：函数值下降曲线
% semilogy(1:maxit, fvalue, 'ob-')
% xlabel('k'); ylabel('|f(x)|');
% title('函数值下降曲线')

%%---------------------------------------------------------------------------------------------------------------

% steffensen加速迭代
% syms x;
% f(x)=(x+1)^(1/3);
% g(x)=x-(f(x)-x)^2/(f(f(x))-2*f(x)+x);
% % in [1,2]
% if(g(1)<1||g(2)>2)
%     fprintf("不收敛");
%     return;
% end
% x=1;
% step=10;
% eps=10^-10;
% while(abs(g(x)-x)>eps&&step>0)
%     step=step-1;
%     fprintf('x=%.5f g=%.5f\n',x,g(x));
%     x=g(x);
% end
% 
% if(abs(g(x)-x)<=eps)
%     fprintf('结果是%.8f\n',x1);
% end
% if(step==0)
%     fprintf('步数不足\n');
% end

%%---------------------------------------------------------------------------------------------------------------

%牛顿法解方程
% syms x;
% f(x)=x^3-x-1;
% df(x)=diff(f,1);
% fprintf('%s\n',char(df));
% n=10;
% eps=10^-15;
% x=0.6;
% for i=1:n
%     x=x-f(x)/df(x);
%     fprintf('k=%d    x=%.6f \n',i,x);
%     if(abs(f(x)/df(x))<eps)
%         break;
%     end
% end

%%---------------------------------------------------------------------------------------------------------------

% 牛顿下山
% syms x;
% f(x)=x^3-x-1;
% df(x)=diff(f,1);
% fprintf('%s\n',char(df));
% n=10;
% eps=10^-15;
% x=0.6;
% r=1;
% while(abs(f(x-r*f(x)/df(x)))>abs(f(x)))
%     r=r/2;
%     x=x-r*f(x)/df(x);
% end
% for i=1:n
%     x=x-f(x)/df(x);
%     fprintf('k=%d    x=%.6f \n',i,x);
%     if(abs(f(x)/df(x))<eps)
%         break;
%     end
% end

%%---------------------------------------------------------------------------------------------------------------

%知道根的重数的牛顿
% syms x;
% f(x)=x^4-4*x^2+4;
% df(x)=diff(f,1);
% step=10;
% eps=10^-7;
% x=1.7;
% x1=x-2*f(x)/df(x);
% for i=1:step
%     fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x1, f(x1));
%     x=x1;
%     x1=x-2*f(x)/df(x);
%     if abs(x1-x)<eps, break, end
% end

%%---------------------------------------------------------------------------------------------------------------

%不知道根的重数的
% syms x;
% f(x)=(x^2-2)^4;
% df(x)=diff(f,1);
% d2f(x)=diff(df,1);
% step=10;
% eps=10^-7;
% x=1.7;
% x1=x-f(x)*df(x)/(df(x)^2-d2f(x)*f(x));
% for i=1:step
%     fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x1, f(x1));
%     x=x1;
%     x1=x-f(x)*df(x)/(df(x)^2-d2f(x)*f(x));
%     if abs(x1-x)<eps, break, end
% end

%%---------------------------------------------------------------------------------------------------------------

%弦截法解方程
% syms x;
% f(x)=x^2-5;
% step=10;
% eps=10^-7;
% x=2.6;
% x1=2.5;
% x2=x1-f(x1)/((f(x1)-f(x))/(x1-x));
% for i=1:step
%     fprintf('k=%d, x=%.8f, f(x)=%.2e\n', i, x2, f(x2));
%     x=x1;
%     x1=x2;
%     x2=x1-f(x1)/((f(x1)-f(x))/(x1-x));
%     if abs(x2-x1)<eps, break, end
% end


%%---------------------------------------------------------------------------------------------------------------

%抛物线法
% syms x;
% f(x) = x^3-x-1;
% x1 = 1;
% x2 = 1.5;
% x3 = 2;
% while 1 %开始迭代
%     y1 = f(x1);
%     y2 = f(x2);
%     y3 = f(x3);
%     f1 = lagrange([x1 x2 x3],[y1 y2 y3]); % 拉格朗日插值出抛物线
%     a = f1(1);b = f1(2); c = f1(3);  % 求出 抛物线的系数 a为二次线系数
% 
%     x4_1 = (-b + sqrt(b^2-4*a*c))/(2*a); % 抛物线求根
%     x4_2 = (-b-sqrt(b^2-4*a*c))/(2*a);
% 
%     if(abs(x1-x4_1)<abs(x1-x4_2))
%         x4 = x4_1;
%     else
%         x4 = x4_2;
%     end
%     if abs(f(x4))< 1e-6
%         if flag ==1||flag==0
%             result = x4;
%         elseif flag ==2
%             result = x4;
%             tol = abs(f(x4));
%         end
%     end
%     x1 = x2;
%     x2 = x3;
%     x3 = x4;
% end


%%---------------------------------------------------------------------------------------------------------------

%高斯消去
% A=[2,3,4;3,5,2;4,3,30];
% b=[6;5;32];
% n=size(A,1);
% for i=1:n-1
%     for j=i+1:n
%         t=A(j,i)/A(i,i);
%         for k=i+1:n
%             A(j,k)=A(j,k)-t*A(i,k);
%         end
%         b(j)=b(j)-t*b(i);
%     end
% end
% x(n)=b(n)/A(n,n);
% for i=n-1:-1:1
%     s=0;
%     for j=i+1:n
%         s=s+A(i,j)*x(j);
%     end
%     x(i)=(b(i)-s)/A(i,i);
% end
% disp(x);

%%---------------------------------------------------------------------------------------------------------------

%现成函数LU
% A=[2,3,4;3,5,2;4,3,30];
% b=[6;5;32];
% % [L,U]=lu(A);
% % y=L\b;
% % x=U\y;
% % disp(x);
% 
%%---------------------------------------------------------------------------------------------------------------

%LU分解
% A=[2,3,4;3,5,2;4,3,30];
% b=[6;5;32];
% n=size(A,1);
% L=eye(n);
% U=zeros(n);
% for i=1:n
%     for j=1:n
%         U(i,j)=A(i,j)-sum(L(i,1:i-1)*U(1:i-1,j));
%         for k=i+1:n
%             L(k,i)=(A(k,i)-sum(L(k,1:i-1)*U(1:i-1,i)))/U(i,i);
%         end
%     end
% end
% 
% x=b;
% for i=2:n
%         x(i)=x(i)-sum(L(i,1:i-1)*x(1:i-1));
% end
% x(n)=x(n)/U(n,n);
% for i=n-1:-1:1
%     x(i)=(x(i)-sum(U(i,i+1:n)*x(i+1:n)))/U(i,i);
% end
% disp(x);

%%---------------------------------------------------------------------------------------------------------------

% 平方根
% %A=[4,2,-2;2,2,-3;-2,-3,14];
% %A=[1,0,0;0,4,0;0,0,9];
% %b=[10;5;4];
% A=[1,1,1,1,1;
%     1,2,2,2,2;
%     1,2,3,3,3;
%     1,2,3,4,4;
%     1,2,3,4,5];
% b=[5;9;12;14;15];
% n=size(A,1);
% for i=1:n
%     A(i,i)=sqrt(A(i,i));
%     A(i+1:n,i)=A(i+1:n,i)/A(i,i);
%     for j=i+1:n
%         A(j:n,j)=A(j:n,j)-A(j:n,i)*A(j,i);
%     end
% end
% 
% %disp(A);
% A=tril(A);% 取下三角A
% for i=1:n-1
%     b(i)=b(i)/A(i,i);
%     b(i+1:n)=b(i+1:n)-b(i)*A(i+1:n,i);
% end
% %disp(b);
% b(n)=b(n)/A(n,n);
% A=A';
% for i=n:-1:2
%     b(i)=b(i)/A(i,i);
%     b(1:i-1)=b(1:i-1)-b(i)*A(1:i-1,i);
% end
% b(1)=b(1)/A(1,1);
% 
%  disp(b);

%%---------------------------------------------------------------------------------------------------------------

%改进平方根
% A=[4,2,-2;2,2,-3;-2,-3,14];
% b=[10;5;4];
% n=size(A,1);
% D=diag(diag(A));
% 
% for i=1:n
%     for j=1:i-1
%         for k=1:j-1
%             A(i,j)=A(i,j)-A(i,k)*A(j,k)*D(k,k);
%         end
%         A(i,j)=A(i,j)/D(j,j);
%     end
%     s=0;
%     for j=1:i-1
%         s=s+A(i,j)*A(i,j)*D(j,j);
%     end
%     D(i,i)=A(i,i)-s;
% end
% %disp(A); %ok
% %disp(D);
% A=A-diag(diag(A))+eye(n);
% A=tril(A);
% %disp(A);
% for i=1:n-1
%     b(i)=b(i)/A(i,i);
%     b(i+1:n)=b(i+1:n)-b(i)*A(i+1:n,i);
% end
% %disp(b);
% b(n)=b(n)/A(n,n);
% A=D*A';
% for i=n:-1:2
%     b(i)=b(i)/A(i,i);
%     b(1:i-1)=b(1:i-1)-b(i)*A(1:i-1,i);
% end
% b(1)=b(1)/A(1,1);
% 
% disp(b);

%%---------------------------------------------------------------------------------------------------------------

%追赶法
% A=[3,1,0;1,2,1;0,1,3];
% b=[3;4;3];
% for i = 2:n
%     A(i,i-1) = A(i,i-1)/A(i-1,i-1);
%     A(i,i) = A(i,i) - A(i-1,i) * A(i,i-1);
%     b(i) = b(i) - b(i-1) * A(i,i-1);
% end
% x(n) = b(n) / A(n,n); 
% for i = n-1 :-1:1
%     x(i) = (b(i) - A(i,i+1) * x(i+1)) / A(i,i);
% end
% 
% disp(x);

%%---------------------------------------------------------------------------------------------------------------

%雅各比迭代法
% A=[9,-1,-1;-1,8,0;-1,0,9];
% b=[7;7;8];
% n=size(A,1);
% x=zeros(n,1);
% eps=1e-7;
% steps=100;
% D=diag(diag(A));
% L=D-tril(A);
% U=D-triu(A);
% 
% for i=1:steps
%     x=D\((L+U)*x+b);
%     rem=norm(b-A*x)/norm(b);
%     fprintf('i=%d, relres=%.2e\n',i, rem);
%     disp(x);
%     if(rem<eps)
%         break;
%     end
% end

%%---------------------------------------------------------------------------------------------------------------

%GS高斯赛德尔迭代
% A=[9,-1,-1;-1,8,0;-1,0,9];
% b=[7;7;8];
% n=size(A,1);
% x=zeros(n,1);
% eps=1e-7;
% steps=100;
% D=diag(diag(A));
% L=D-tril(A);
% U=D-triu(A);
% 
% for i=1:steps
%     x=(D-L)\(U*x+b);
%     rem=norm(b-A*x)/norm(b);
%     fprintf('i=%d, relres=%.2e\n',i, rem);
%     disp(x);
%     if(rem<eps)
%         break;
%     end
% end

%%---------------------------------------------------------------------------------------------------------------

% SOR 超松弛迭代
% A=[9,-1,-1;-1,8,0;-1,0,9];
% b=[7;7;8];
% n=size(A,1);
% x=zeros(n,1);
% eps=1e-7;
% steps=100;
% D=diag(diag(A));
% L=D-tril(A);
% U=D-triu(A);
% 
% w=1;
% for i=1:steps
%     x=(D-w*L)\((w*U+(1-w)*D)*x+w*b);
%     rem=norm(b-A*x)/norm(b);
%     fprintf('i=%d, relres=%.2e\n',i, rem);
%     disp(x);
%     if(rem<eps)
%         break;
%     end
% end

%%---------------------------------------------------------------------------------------------------------------

%拟合
% x1 = [1 2 3 4 5];
% y1 = [4  4.5  6  8  8.5];
% subplot(1,3,1);
% p=polyfit(x1,y1,1);
% X1=0:0.01:6;
% Y1=p(1)*X1+p(2);
% plot(x1,y1,'o',X1,Y1,'r-');
% xlabel('x');ylabel('y');
% legend('data','y=ax+b','Location','northwest');
%
% x2 = [1.00  1.25  1.50  1.75  2.00];
% y2 = [5.10  5.79  6.53  7.45  8.46];
% subplot(1,3,2);
% y21=log(y2);
% p=polyfit(x2,y21,1);
% X2=0:0.01:4;
% Y2=exp(p(1)*X2+p(2));
% plot(x2,y2,'o',X2,Y2,'r-');
% xlabel('x');ylabel('y');
% legend('data','y=e^(ax+b)','Location','northwest');
% 
% x3 = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
% y3 = [1, 1.75, 1.96, 2.19, 2.44, 2.71, 3.00];
% subplot(1,3,3);
% p=polyfit(x3,y3,2);
% X3=0:0.01:2;
% Y3=polyval(p,X3);
% plot(x3,y3,'o',X3,Y3,'r-');
% xlabel('x');ylabel('y');
% legend('data','y=ax^2+bx+c','Location','northwest');

%%---------------------------------------------------------------------------------------------------------------

%Hermite 三点三次埃尔米特插值
% syms x;
% x0=1/4;
% x1=1;
% x2=9/4;
% f(x)=(x)^(3/2);
% df(x)=diff(f,x);
% a01=(f(x0)-f(x1))/(x0-x1);
% a12=(f(x2)-f(x1))/(x2-x1);
% a20=(a12-a01)/(x2-x0);
% alpha=(df(x1)-a01-a20*(x1-x0))/((x1-x0)*(x1-x2));
% fprintf("%.6f",alpha);
% y(x)=f(x0)+a01*(x-x0)+a20*(x-x0)*(x-x1)+alpha*(x-x0)*(x-x1)*(x-x2);
% fprintf("%s",char(y));

%%---------------------------------------------------------------------------------------------------------------

%龙格现象
% f = @(x) 1./(1+x.^2);
% a = -5; 
% b = 5;
% 
% xi = -5: 0.01 :5;  % 用于绘图的等距离散点
% yt = f(xi);  % 函数在这些点上的精确值
% 
% for n = 2  : 2 : 12
%     % n = 5;
%     X = a : (b-a)/n : b;  % 插值节点
%     Y = f(X);
% 
%     yi = zeros(1,length(xi));  % 存储通过插值计算的近似值
%     for i = 1 : length(xi)
%         for k = 0 : n
% 	    xtmp1 = xi(i)-X([1:k,k+2:end]);
% 	    xtmp2 = X(k+1)-X([1:k,k+2:end]);
%             yi(i) = yi(i) + Y(k+1)*prod(xtmp1./xtmp2);
%         end
%     end
% 
%     subplot(1,6,n/2);
% 
%     plot(xi,yt,'r-', xi,yi,'b-','LineWidth',2);
%     hold on;
%     plot(X,Y,'bo','LineWidth',2,'markersize',10);
%     hold off;
%     axis([-5,5,-4,2]);
%     tit = ['n=',int2str(n)];
%     title(tit,'FontSize',20); 
%     legend({'f(x)','L_n(x)'},'FontSize',20,'FontName', 'Times New Roman','FontAngle','italic');
% end

%%---------------------------------------------------------------------------------------------------------------

%Hermite 两点三次埃尔米特插值（分段三次埃尔米特插值，分段低次）
% f=@(x)1./(1+x.^2);
% x=-5:0.01:5;
% y=f(x);
% df = @(x) -(2*x)./(x.^2 + 1).^2; 
% plot(x,y);
% hold on;
% for i=-5:1:4
%     x0=i;
%     x1=i+1;  
%     y1=@(x)...
%             f(x0).*(1+2.*(x-x0)./(x1-x0)).*((x-x1)./(x0-x1)).^2+...
%             f(x1).*(1+2.*(x-x1)./(x0-x1)).*((x-x0)./(x1-x0)).^2+...
%             df(x0).*(x-x0).*((x-x1)./(x0-x1)).^2+...
%             df(x1).*(x-x1).*((x-x0)./(x1-x0)).^2;
%     N=500;
%     x2=x0:(x1-x0)/N:x1;
%     y2=y1(x2);  
%     plot(x2,y2,'r-');
% end
% hold off;
%--------------------------------
% clear; clc; 
% close all;
% 
% f = @(x) 1./(1+x.^2); % 函数表达式
% df = @(x) -(2*x)./(x.^2 + 1).^2; % 一阶导数, 用于Hermite插值
% a = -5;  b = 5;  % 插值区间
% n = 10;          % 区间等分数
% h = (b-a)/n;     % 步长 
% xi = a : h : b;  % 插值节点
% fi = f(xi);      % 插值节点上的函数值
% dfi = df(xi);    % 一阶导数值
% 
% x = a : (b-a)/500 : b; % 需要插值的点，用于绘图
% 
% % 定义线性插值函数
% L1 = @(x,x0,x1,f0,f1) f0*(x-x1)/(x0-x1) + f1*(x-x0)/(x1-x0);
% 
% % 定义两点三次Hermite插值函数
% H3 = @(x,x0,x1,f0,f1,df0,df1) ...
%    (f0*(1+2*(x-x0)/(x1-x0))+df0*(x-x0))*((x-x1)/(x0-x1))^2 + ...
%    (f1*(1+2*(x-x1)/(x0-x1))+df1*(x-x1))*((x-x0)/(x1-x0))^2;
% 
% % 分段插值
% N = length(x);
% y1 = zeros(1,N);  % 分段线性插值
% y2 = zeros(1,N);  % 分段三次Hermite插值
% for j = 1 : N
%    for k=1:n+1   % 寻找 x(j) 所在的小区间 [x_k,x_{k+1}]
%       if xi(k) >= x(j)
%          break;  % 找到第一个不小于 x(j) 的插值节点
%       end
%    end
%    if k>1
%       k=k-1;
%    end
% %    fprintf('x(j)=%.2f, k=%d, 插值小区间: [%.2f,%.2f]\n', ...
% %       x(j), k, xi(k),xi(k+1))
%    y1(j) = L1(x(j),xi(k),xi(k+1),fi(k),fi(k+1));
%    y2(j) = H3(x(j),xi(k),xi(k+1),fi(k),fi(k+1),dfi(k),dfi(k+1));
% end
% 
% % 绘图
% hold on;
% plot(x,f(x),'r-', 'LineWidth',2);  % f(x) 图像
% plot(x,y1,'b-', 'LineWidth',2,'markersize',5); % 分段线性插值图像
% plot(x,y2,'k-', 'LineWidth',2,'markersize',10); % 分段三次Hermite插值图像
% titstr=['n=',int2str(n)]; title(titstr,'fontsize',14);
% axis([-5,5,-1,2]);
% 
% plot(xi,fi,'ok','markersize',10,'LineWidth',1.5); % 绘制插值节点
% legend('f(x)','piecewise L1(x)','piecewise H3(x)','data');
% set(gca,'FontSize',15);
% hold off

%%---------------------------------------------------------------------------------------------------------------

% 分段线性插值
% f=@(x)1./(1+x.^2);
% X=-5:0.01:5;
% Y=f(X);
% plot(X,Y);
% hold on;
% for i=-5:4
%     x0=i;
%     x1=i+1;
%     x2=x0:0.01:x1;
%     y2=f(x0)*(x2-x1)/(x0-x1)+f(x1)*(x2-x0)/(x1-x0);
%     plot(x2,y2,'r-');
% end
% hold off;
