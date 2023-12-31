%
% Example 1：一次多项式拟合
%
clc; clear; close all;
x = [1 2 3 4 5];
y = [4  4.5  6  8  8.5];
n = 1;
p = polyfit(x,y,n);
xx = [0:0.01:6];
S = p(1)*xx + p(2);
plot(x,y,'o',xx,S,'-')
xlabel('x'); ylabel('y')
legend('data','y=ax+b','Location','northwest')
fprintf('%d 次最小二乘拟合多项式为: %s\n',n, char(poly2sym(p)));
%
% Example 2：指数函数拟合
%
clear;
x = [1.00  1.25  1.50  1.75  2.00];
y = [5.10  5.79  6.53  7.45  8.46];
y1=log(y);
% 1次拟合；
p=polyfit(x,y1,1);
xx=[0.5:0.01:2.5];
S=exp(p(2))*exp(p(1)*xx);
figure;
plot(x,y,'o',xx,S,'-')
xlabel('x'); ylabel('y')
legend('data','y=ae^{bx}','Location','northwest')
%
% Example 3：二次多项式拟合
%
clear;
x = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
y = [1, 1.75, 1.96, 2.19, 2.44, 2.71, 3.00];
n = 2;
p=polyfit(x,y,n);
fprintf('%d 次最小二乘拟合多项式为: %s\n',n, char(poly2sym(p)));
% 绘图
figure
xx = 0:0.01:1;
plot(x, y, 'o',xx, polyval(p,xx),'-');
xlabel('x'); ylabel('y')
%plot(x,polyval(p,x),'bo','MarkerSize',10); % 绘拟合曲线在给定数据点上的近似值