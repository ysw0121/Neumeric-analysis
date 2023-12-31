% clear all;
% close all;

A=input('输入线性方程组矩阵A=');
b=input('输入结果列向量b, AX=b, b=');
n=size(A,1);
x0=zeros(n,1);
eps=1e-5;
maxTimes=100;

D = diag(diag(A));  % A 的对角线部分
L = D - tril(A);    % -L 为 A 的严格下三角部分
U = D - triu(A);    % -U 为 A 的严格上三角部分

%Jacobi
fprintf('Jacobi: \n')
x1 = x0;  
for i= 1 : maxTimes
    x1 = D \ ( (L+U)*x1 + b );
    relres = norm(b - A*x1) / norm(b); % norm--范数
    fprintf('i=%d, relres=%.2e',i, relres);
    fprintf('x=');
    for j=1:n
        fprintf('%.4f  ',x1(j));
    end
    fprintf('\n');
    if (relres<eps), break, end
end

% G-S
fprintf('\nG-S: \n')
x2 = x0; 
for i = 1 :maxTimes
    x2 = (D-L) \ ( U*x2 + b );
    relres = norm(b - A*x2) / norm(b); 
    fprintf('i=%d, relres=%.2e',i, relres)
    fprintf('x=');
    for j=1:n
        fprintf('%.4f  ',x2(j));
    end
    fprintf('\n');
    if (relres<eps), break, end
end

% SOR
fprintf('\nSOR: \n')
omega = 1.3;  
x3 = x0; 
for i= 1 : maxTimes
    x3 = (D-omega*L) \ ( ((1-omega)*D + omega*U)*x3 + omega*b );
    relres = norm(b - A*x3) / norm(b); % 相对残量
    fprintf('i=%d, relres=%.2e',i, relres)
    fprintf('x=');
    for j=1:n
        fprintf('%.4f  ',x3(j));
    end
    fprintf('\n');
    if (relres<eps), break, end
end

