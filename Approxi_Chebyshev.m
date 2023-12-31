% 对Chebyshev多项式画图
%
clear
% 多项式次数从0到N
N = 6;
syms T x
% 画图曲线颜色设置
c={'b','r','g','k','m','c','y'};
for n = 0:N
    % Tn表达式
    T = cos(n*acos(x));
    h = fplot(T,[-1 1],'Linewidth',1);
    set(h,'Color',c{n+1})
    hold on
end
legend('n=0','n=1','n=2','n=3','n=4','n=5','n=6','Location','NorthEastOutside')
title('Chebyshev polynomials');
xlabel('x')