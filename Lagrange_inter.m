% homework 1-------Lagrange inter for expected lifetime


%   output:   the lifetime at 70℃ is  49.7×1000 hrs

x1=25; x2=40; x3=50; x4=60;
y1=95; y2=75; y3=63; y4=54;

x=70;

y=y1*(x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))+y2*(x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))+y3*(x-x2)*(x-x1)*(x-x4)/((x3-x2)*(x3-x1)*(x3-x4))+y4*(x-x2)*(x-x3)*(x-x1)/((x4-x2)*(x4-x3)*(x4-x1));

fprintf("the lifetime at 70℃ is  %.1f×1000 hrs",y);