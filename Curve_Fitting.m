%
% Example 1��һ�ζ���ʽ���
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
fprintf('%d ����С������϶���ʽΪ: %s\n',n, char(poly2sym(p)));
%
% Example 2��ָ���������
%
clear;
x = [1.00  1.25  1.50  1.75  2.00];
y = [5.10  5.79  6.53  7.45  8.46];
y1=log(y);
% 1����ϣ�
p=polyfit(x,y1,1);
xx=[0.5:0.01:2.5];
S=exp(p(2))*exp(p(1)*xx);
figure;
plot(x,y,'o',xx,S,'-')
xlabel('x'); ylabel('y')
legend('data','y=ae^{bx}','Location','northwest')
%
% Example 3�����ζ���ʽ���
%
clear;
x = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
y = [1, 1.75, 1.96, 2.19, 2.44, 2.71, 3.00];
n = 2;
p=polyfit(x,y,n);
fprintf('%d ����С������϶���ʽΪ: %s\n',n, char(poly2sym(p)));
% ��ͼ
figure
xx = 0:0.01:1;
plot(x, y, 'o',xx, polyval(p,xx),'-');
xlabel('x'); ylabel('y')
%plot(x,polyval(p,x),'bo','MarkerSize',10); % ����������ڸ������ݵ��ϵĽ���ֵ