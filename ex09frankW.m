%ʹ��Frank-Wolfe�����
clear;
clc;

f=@(x1,x2) -20*x1-16*x2+2*x1^2+x2^2+(x1+x2)^2;
x0=[0; 0];

minx=frankWolfe(f, x0, 0.25);

% x=-5: 0.01: 15; %����������άͼ�β鿴һ�¿������ڴ�ŵ���״
% y=-5: 0.01: 10;
% [X,Y]=meshgrid(x,y);
% Z=-(20*X+16*Y-2*X.^2-Y.^2-(X+Y).^2);
% % plot3(X,Y,Z);
% contour3(X,Y,Z,80); %���ȸ���ͼ