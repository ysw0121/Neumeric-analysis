%杜利特尔分解法解方程

% clear all;
% close all;


A=input('输入线性方程组矩阵A=');
b=input('输入结果列向量b, AX=b, b=');
n=size(A,1);

for i=1:n
    pd=A(1:i,1:i);
    if det(pd)==0
        disp('Doolittle分解不存在'); %前n-1阶顺序主子式均不能为0
        return;
    end
end

L=eye(n);
U=zeros(n);
for i=1:n
    for j=i:n
        U(i,j)=A(i,j)-sum(L(i,1:i-1)*U(1:i-1,j));
        for k=i+1:n
            L(k,i)=(A(k,i)-sum(L(k,1:i-1)*U(1:i-1,i)))/U(i,i);
        end
    end
end


% disp(L);
% disp(U);


x=b;
% disp(x);


%Ly=b
for i=2:n
    for j=1:i-1
    x(i)=x(i)-L(i,j)*x(j);
    end
end
%disp(x);

%Ux=y
x(n)=x(n)/U(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+U(i,j)*x(j);
    end
    x(i)=(x(i)-s)/U(i,i);
end

fprintf('det(A)=%.4f\n',det(U));
disp('X=');
disp(x);