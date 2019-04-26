function minx=frankWolfe(f, x0, eps)
% minx�����������Ľ��
% f����Ҫ���ĺ���
% x0����ʼֵ
% eps����⾫��

syms x1 x2 m;
d=[diff(f,x1); diff(f,x2)]; %�ݶȷ���

dx0=subs(d, x1, x0(1));     %��ֵ���������ʼֵ������ݶ�ֵ
dx0=double( subs(dx0, x2, x0(2))); 

fy=dx0' ;
A=[1 1; -1 0; 0 -1];
b=[5 0 0];
y0=linprog(fy, A, b);

tol=abs( dx0' * (y0-x0));

while tol>eps
    xtmp=x0+m*(y0-x0);  %����ֵ�����һԪ��Сֵ����
    fm=subs(sym(f), x1, xtmp(1));
    fm=subs(fm, x2, xtmp(2));  %fm ����m��һԪ����
    
    m0=newtonMin(fm, 0, eps); %������ʣ�ǰ������
    
    x01=x0+m0*(y0-x0); %��ǰ���������Ĺ�ʽ
    
    x0=x01;
    
    dx0=subs(d, x1, x0(1));     %��ֵ���������ʼֵ������ݶ�ֵ
    dx0=double( subs(dx0, x2, x0(2)));
    
    fy=dx0' ;
    A=[1 1; -1 0; 0 -1];
    b=[5 0 0];
    y0=linprog(fy, A, b);

    tol=abs( dx0' * (y0-x0));
end
minx=x0;

function minx=newtonMin(f,x0,eps)
% ţ�ٷ�����f�ϵ���Сֵ
% f:������
% x0=7; %��x0��ʼ����
% eps����⾫��
% minx:��

df=diff(sym(f));    %����
ddf=diff(sym(df));  %���ε���

dfx0=subs(sym(df),symvar(sym(df)),x0);
ddfx0=subs(sym(ddf),symvar(sym(ddf)),x0);

x1=double ( x0-dfx0/ddfx0 );  %ţ�ٷ������ĵ�����ʽ
tol=double(abs(x1-x0));

while tol>eps
    x0=x1;
    dfx0=subs(sym(df),symvar(sym(df)),x0);
    ddfx0=subs(sym(ddf),symvar(sym(ddf)),x0);
    x1=double(x0-dfx0/ddfx0);
    tol=double(abs(x1-x0));
end
minx=x1;

