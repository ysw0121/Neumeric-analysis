%   problem(4)

init=100;
answer=0;
for i=1:10
    answer=answer+2*init;
    init=init*0.5;
end
answer=answer-100;
disp(answer);
disp(init);

% output:   sum= 299.6094,  the 10th height = 0.0977