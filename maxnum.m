%  problem(2)
a=input('pls input 4*5 matrix a= ');
[x,y]=find(a==max(max(a)));
disp(max(max(a)));
disp([x,y]);

% e.g.   input:  [1,2,2,3,34;5,5,6,7,45;9,5,7,8,1;12,43,23,66,44]
%        output:   66
%                   4     4