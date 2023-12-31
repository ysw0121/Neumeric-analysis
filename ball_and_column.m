r=2;
dim=4;
t=linspace(-2,2);
q=linspace(-2,2);

[x,y,z]=meshgrid(t,q/2,q);
c=cat(dim,x.^2+y.^2+z.^2-r^2,x.^2+y.^2-r*x);
% disp(c);
v=max(c,[],dim);
% disp(v);

subplot(3,3,1);
isosurface(x,y,z,v,0); % 取等值面
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(20,15);
lighting phong;  
grid on;

subplot(3,3,2);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(20,-5);
lighting phong;  
grid on;

subplot(3,3,3);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(110,10);
lighting phong;  
grid on;

subplot(3,3,4);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(110,-10);
lighting phong;  
grid on;

subplot(3,3,5);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(190,15);
lighting phong;  
grid on;

subplot(3,3,6);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(190,-15);
lighting phong;  
grid on;

subplot(3,3,7);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(330,20);
lighting phong;  
grid on;

subplot(3,3,8);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(330,0);
lighting phong;  
grid on;

subplot(3,3,9);
isosurface(x,y,z,v,0);
% colormap cool;
colormap([0.5,0.9,0]);   
brighten(0.5);       
camlight right;
view(120,0);
lighting phong;  
grid on;