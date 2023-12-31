%cardioid



subplot(2,2,1);
theta=linspace(0,2*pi);
r=3*(1-cos(theta));  % ρ=a(1-cosθ)
polarplot(theta,r,'m','LineWidth',1);
title('心形线：ρ=a(1-cosθ)，a=3');
% gtext('ρ=a(1-cosθ)');
legend('line1');
grid on;

%Hyperbolic paraboloid
subplot(2,2,2);
x=-5:0.1:5;
y=-5:0.1:5;
Axis=([-5,5,-5,5]);
[X,Y]=meshgrid(x,y);
z=X.^2/3-Y.^2/3;
mesh(X,Y,z);
title('马鞍面z=x^2/a-y^2/b，a,b=3');
xlabel('variable x');
ylabel('variable y');
zlabel('variable z');
legend('line2');
%shading flat;
colormap cool;
%material shiny;
grid on;


subplot(2,2,3);
mesh(X,Y,z);
view([110,20]);
title('马鞍面视角1');
xlabel('variable x');
ylabel('variable y');
zlabel('variable z');
legend('line3');


subplot(2,2,4);
mesh(X,Y,z);
view([175,20]);
title('马鞍面视角2');
xlabel('variable x');
ylabel('variable y');
zlabel('variable z');
legend('line4');