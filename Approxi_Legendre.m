% 对Legendre多项式画图
% 
clear
% 多项式次数从0到N
N = 6;
syms P x
% 画图曲线颜色设置
c={'b','r','g','k','m','y','c'};
for n = 0:N
    % Legendre多项式的罗德里克表达式（也可用递推关系计算）
    P = 1/(2^n*prod(1:n))*diff((x^2-1)^n,x,n);
    h = fplot(P,[-1 1],'LineWidth',1);
    set(h,'Color',c{n+1})
    hold on
end
legend('n=0','n=1','n=2','n=3','n=4','n=5','n=6','Location','NorthEastOutside')
title('Legendre polynomials');
xlabel('x')