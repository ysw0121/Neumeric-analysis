%clear all,clc;
% t=-0.5:.01:0.5;
% [x,y]=meshgrid(t);%形成格点矩阵
% z=sin(4*pi*x)+cos(6*pi*y);
% figure(1)
% mesh(x,y,z);
% axis([-0.5 0.5 -0.5 0.5 -4 4]);
% title('z')
% colormap cool%cool是一种配色方案，还有其他方案如winter，summer····见help colormap
% colorbar


% xt = @(t) cos(3*t);
% yt = @(t) sin(2*t);
% fplot(xt,yt)
% fplot (@(x) tanh(x),[-2,2])

% fplot(@(x) exp(x),[-3 0],'b')
% hold on                %原图保留，在此基础上做函数图像
% fplot(@(x) cos(x),[0 3],'r')
% hold off
% grid on  %开网格


% fplot(@(x) sin(x+pi/5),'Linewidth',2);
% hold on
% fplot(@(x) sin(x-pi/5),'--or');
% fplot(@(x) sin(x),'-.*c')
% hold off

% x=0:.1:10;
% semilogy(x,10.^x)

% x=logspace(-1,2);
% loglog(x,exp(x),'-s')
% grid on

% x=[1:1:100];
% subplot(2,3,1);
% plot(x,x.^3);
% grid on;
% title 'plot-y=x^3';
% 
% subplot(2,3,2);
% loglog(x,x.^3);
% grid on;
% title 'loglog-logy=3logx';
% 
% subplot(2,3,3);
% plotyy(x,x.^3,x,x);
% grid on;
% title 'plotyy-y=x^3,logy=3logx';
% 
% 
% 
% subplot(2,3,4);
% semilogx(x,x.^3);
% grid on;
% title 'semilogx-y=3logx';
% 
% subplot(2,3,5);
% semilogy(x,x.^3);
% grid on;
% title 'semilogy-logy=x^3';
% hold off




% x=-3:0.1:3;
% y=1:0.1:5;
% [X,Y]=meshgrid(x,y);
% Z=(X+Y).^2;
% surf(X,Y,Z)
% shading flat

% x=-3:0.1:3;   y=1:0.1:5;
%      [X,Y]=meshgrid(x,y);
%      Z=(X+Y).^2;
%      mesh(X,Y,Z)     


% x=-3:0.1:3;   y=1:0.1:5;
%      [X,Y]=meshgrid(x,y);
%      Z=(X+Y).^2;
%      subplot(2,2,1); mesh(X,Y,Z)
%      subplot(2,2,2);mesh(X,Y,Z);view(50,-34)
%      subplot(2,2,3);mesh(X,Y,Z);view(-60,70) 
%      subplot(2,2,4);mesh(X,Y,Z);view([0,1,1])


% [x,y,z]=peaks(30);
% surf(x,y,z)
% axis([-3 3 -3 3 -10 10])
% %axis off
% %shading interp
% %colormap(hot)
% m=moviein(15);
% for i=1:15
%    view(-37.5+24*(i-1),30)
%    m(:,i)=getframe;
% end
% movie(m)

% [x,y,z]=peaks;
%         subplot(1,2,1)       
%         contour3(x,y,z,16,'s')   
%         grid,   xlabel('x-axis'),  ylabel('y-axis')
%         zlabel('z-axis')
%         title('contour3 of peaks'); 
%         subplot(1,2,2)
%         contour(x,y,z,16,'s')
%         grid,   xlabel('x-axis'),    ylabel('y-axis')
%         title('contour of peaks');
% 
% [x,y,z]=sphere(16);
% X=[x(:)*.5 x(:)*.75 x(:)];
% Y=[y(:)*.5 y(:)*.75 y(:)];
% Z=[z(:)*.5 z(:)*.75 z(:)];
% S=repmat([1 .75 .5]*10,prod(size(x)),1);
% C=repmat([1 2 3],prod(size(x)),1);
% scatter3(X(:),Y(:),Z(:),S(:),C(:),'filled'),view(-60,60)



% [X,Y]=meshgrid(-2:.2:2,-2:.2:3);
% Z=X.*exp(-X.^2-Y.^2);
% subplot(1,2,1);
% [C,h]=contour(X,Y,Z);
% clabel(C,h)
% colormap cool
% subplot(1,2,2);
% contour3(X,Y,Z,'m');


theta=linspace(0,2*pi),                         
rho=sin(2*theta).*cos(2*theta);
polar(theta,rho,'g');
title('Polar plot of sin(2*theta).*cos(2*theta)');