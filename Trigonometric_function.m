clc; clear; close all;
% cos x与sin x图
x=linspace(-pi,pi,50); 
y1=cos(x);
y2=sin(x);
plot(x,y1,x,y2,'LineWidth',2);
xlim([-1.5*pi,1.5*pi]);
ylim([-1.5,1.5]);
set(gca,'xtick',-pi:pi/2:pi);%设置x轴坐标间隔
set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
ax1 = gca;
ax1.XAxisLocation = 'origin';
ax1.YAxisLocation = 'origin';
legend({'cos(x)','sin(x)'},'FontSize',10);

