%problem(3)
function [answer]=factsum(n)
answer=0;
for i=1:n
    fac=1;
    for j=1:i
        fac=fac*j;
    end
    answer=answer+fac;
end
end

%>>factsum(20)
%   ans =
% 2.5613e+18
