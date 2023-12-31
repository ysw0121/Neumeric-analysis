% with Newton

% output: 2010 oil production is -1951646.134×10^6 bbl/day 

% 由结果可知，结果明显不符合常识，而且偏差巨大。
% 因此，结合图像，本插值方法出现Runge现象。

% 此模型中，9次多项式不是一个很好的方法。
% 本题在 x=11，即下一个点的时候就出现了巨大偏差，无法用于预测，
% 而且本次数据分布不算均匀，因此产生误差很大。


% 虽然牛顿插值公式具有高精度的特点，可以在一定程度上减小插值误差。
% 而且计算过程相对简单，容易实现。
% 并且可以通过递推的方式来计算多项式系数。在插值点数较多时，计算效率较高。
% 但是，牛顿插值公式的多项式次数会随着插值点数的增加而增加，
% 因此在插值点数较多时，多项式函数的次数可能会非常高，导致计算复杂度增加。
% 而且牛顿插值公式对于插值点的分布比较敏感，如果插值点分布不均匀，可能会导致插值误差较大。



fprintf('2010 oil production is %.3f×10^6 bbl/day \n',insertion(17));

X=1:10;
Y=[67.052, 68.008, 69.803, 72.024, 73.400, 72.063, 74.669, 74.487, 74.065, 76.777];

subplot(1,4,1);
scatter(X,Y,'filled','b');
hold on;
fplot(@(x) insertion(x),[1,10]);
hold off;

subplot(1,4,2);
scatter(X,Y,'filled','b');
hold on;
fplot(@(x) insertion(x),[1,11]);
hold off;


subplot(1,4,3);
scatter(X,Y,'filled','b');
hold on;
fplot(@(x) insertion(x),[1,14]);
hold off;

subplot(1,4,4);
scatter(X,Y,'filled','b');
hold on;
fplot(@(x) insertion(x),[1,20]);
hold off;



function [p, q] = divided_difference(x,y)
m = length(x);
x = x(:);
p = zeros(m, m+1);
p(:,1) = x; 
p(:,2) = y(:);
for k = 3 : m+1
    p(k-1:m, k) = diff(p(k-2:m, k-1)) ./ ( x(k-1:m) - x(1:m+2-k) );
end
q = diag(p(1:m,2:m+1));
end


function[y0]= insertion(x0)
X=1:10;
Y=[67.052, 68.008, 69.803, 72.024, 73.400, 72.063, 74.669, 74.487, 74.065, 76.777];

[p,q]=divided_difference(X,Y);
y0=q(1);

for i=1:9
    a=q(i+1);
    for j=1:i
        a=a*(x0-X(j));
    end
    y0=y0+a;
end
end