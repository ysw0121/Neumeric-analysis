%   problem(1)

a=input('pls input 10 numbers a=');
for n=1:9
    for m=1:10-n
        if a(m)>a(m+1)
            temp=a(m);
            a(m)=a(m+1);
            a(m+1)=temp;
        end
    end
end
disp(a);

%  e.g. input : [10, 4, 3, 5, 88, 50, 91, 52, 1, 3]
%       output:  1     3     3     4     5    10    50    52    88    91
