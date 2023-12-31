%高斯消元法解线性方程
%

% clear all;
% close all;

A=input('输入线性方程组矩阵A=');
b=input('输入结果列向量b, AX=b, b=');
n=size(A,1);
for i=1:n-1
    for j=i+1:n
        t=A(j,i)/A(i,i);
        for k=i+1:n
            A(j,k)=A(j,k)-t*A(i,k);
        end
        b(j)=b(j)-t*b(i);
    end
end
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+A(i,j)*x(j);
    end
    x(i)=(b(i)-s)/A(i,i);
end
disp('X=')
disp(x);
