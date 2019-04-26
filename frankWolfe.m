function minx=frankWolfe(f, x0, eps)
% minx：函数求解出的结果
% f：需要求解的函数
% x0：初始值
% eps：求解精度

syms x1 x2 m;
d=[diff(f,x1); diff(f,x2)]; %梯度方向

dx0=subs(d, x1, x0(1));     %数值化，带入初始值，求出梯度值
dx0=double( subs(dx0, x2, x0(2))); 

fy=dx0' ;
A=[1 1; -1 0; 0 -1];
b=[5 0 0];
y0=linprog(fy, A, b);

tol=abs( dx0' * (y0-x0));

while tol>eps
    xtmp=x0+m*(y0-x0);  %带入值，变成一元最小值问题
    fm=subs(sym(f), x1, xtmp(1));
    fm=subs(fm, x2, xtmp(2));  %fm 关于m的一元函数
    
    m0=newtonMin(fm, 0, eps); %求出最适，前进长度
    
    x01=x0+m0*(y0-x0); %向前迭代，核心公式
    
    x0=x01;
    
    dx0=subs(d, x1, x0(1));     %数值化，带入初始值，求出梯度值
    dx0=double( subs(dx0, x2, x0(2)));
    
    fy=dx0' ;
    A=[1 1; -1 0; 0 -1];
    b=[5 0 0];
    y0=linprog(fy, A, b);

    tol=abs( dx0' * (y0-x0));
end
minx=x0;

function minx=newtonMin(f,x0,eps)
% 牛顿法求函数f上的最小值
% f:函数名
% x0=7; %从x0开始迭代
% eps：求解精度
% minx:解

df=diff(sym(f));    %导数
ddf=diff(sym(df));  %二次导数

dfx0=subs(sym(df),symvar(sym(df)),x0);
ddfx0=subs(sym(ddf),symvar(sym(ddf)),x0);

x1=double ( x0-dfx0/ddfx0 );  %牛顿法，核心迭代公式
tol=double(abs(x1-x0));

while tol>eps
    x0=x1;
    dfx0=subs(sym(df),symvar(sym(df)),x0);
    ddfx0=subs(sym(ddf),symvar(sym(ddf)),x0);
    x1=double(x0-dfx0/ddfx0);
    tol=double(abs(x1-x0));
end
minx=x1;

