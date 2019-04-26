%使用Frank-Wolfe法求解
clear;
clc;

f=@(x1,x2) -20*x1-16*x2+2*x1^2+x2^2+(x1+x2)^2;
x0=[0; 0];

minx=frankWolfe(f, x0, 0.25);

% x=-5: 0.01: 15; %用来绘制三维图形查看一下可行域内大概的形状
% y=-5: 0.01: 10;
% [X,Y]=meshgrid(x,y);
% Z=-(20*X+16*Y-2*X.^2-Y.^2-(X+Y).^2);
% % plot3(X,Y,Z);
% contour3(X,Y,Z,80); %画等高线图