%FractionalFPK.m
% ������ʹ�������ֵС�����Է�����FPK���̽�����⡣
% 5/5/2019

%------------------------------------------------------------------------
function FractionalFPK() %main function
close all;
clear; clc;

% �����޸ĵ���-----------------------------------------------------------
global W0 W1 W2

j = 6;
r = 3;
interval = 1;      %0:���÷�����С��,1:��������С��

alfa = 0.8; % ������
%------------------------------------------------------------------------
xmin = 0.0; xmax = 1.0;
tmin = 0; tmax = 0.01; tao = 0.0001;


deltax = (xmax - xmin) / 2 ^ j;

i = 0 : 2 ^ j;

X = xmin : deltax : xmax;
U0 = X .* (1 - X); % ��ֵ

Uexact = X .* (1 - X) + (2 * X - 3) * tmax ^ alfa / gamma(1 + alfa) - ...
    2 * tmax ^ (2 * alfa) / gamma(1 + 2 * alfa);

num = 0;

dimnum = 2 ^ j + 1;
if interval == 0  % ���÷�����С��
    [W0,W1,W2]=dxfai0(j,i,r,xmin,xmax,X); %theta��һ�׺Ͷ��׵���ֵ(������С��)
else
    [W0,W1,W2] = interval012d(j,r,xmin,xmax,X);%Ȩ������������(����С��)
end

A = tao ^ alfa * (-W2 + W1) + diag(ones(dimnum,1));
invA = inv(A);

n = 0;
U(n + 1,:) = U0;
g(n + 1) = 1; 

n = n + 1;
len = 2 ^ j + 1;
Ut = reshape(U(n,:),len,1);
U(n + 1,:) = invA * Ut;
% boundary condition
U(n + 1,1) = -3 * tao ^ alfa / gamma(1 + alfa) - 2 * tao ^ (2 * alfa) / gamma(1 + 2 * alfa);
U(n + 1,len) = -tao ^ alfa / gamma(1 + alfa) - 2 * tao ^ (2 * alfa) / gamma(1 + 2 * alfa);
g(n + 1) = -alfa * g(n);

for t = tmin + tao : tao : tmax - tao
    temp = 0;
    for k = 0 : n
        temp = temp + g(k + 1) * U(1,:);
    end
    for k = 1 : n
        temp = temp - g(k + 1) * U(n - k + 2,:);
    end
    
    n = n + 1;
    g(n + 1) = (1 - (1 + alfa) / n) * g(n);
    U(n + 1,:) = invA * reshape(temp,len,1);
    % boundary condition
    U(n + 1,1) = -3 * (t + tao) ^ alfa / gamma(1 + alfa) - 2 * (t + tao) ^ (2 * alfa) / gamma(1 + 2 * alfa);
    U(n + 1,len) = -(t + tao) ^ alfa / gamma(1 + alfa) - 2 * (t + tao) ^ (2 * alfa) / gamma(1 + 2 * alfa);
end
figure;

subplot(2,1,1);
plot(X,Uexact,'k-',...
    X,U(n + 1,:),'b-o',...
'LineWidth',2);

abs_error = abs(U(n + 1,:) - Uexact);
subplot(2,1,2); % �������ͼ
plot(X,abs_error,'k-','LineWidth',2);

max_error = max(abs_error)
error_2 = sqrt(sum((U(n+1,:) - Uexact).^2))
%------------------end of the main function--------------------------------

%=============================����С��=================================
% ��������С������ϵ��
% ����k��m�����������,kȷ����,mȷ����
% valueӦ���Ǿ���
function value = amk(j,k,m,L)
% ������kת���ɸ�����ͬ�ľ���������m��ͬ
t = ones(1,length(m));
if(size(k,1) == 1) %��kת����������
    k = k';
end
k = k * t;
%��������mת���ɸ�����ͬ�ľ�������ͬk����
t = ones(size(k,1),1);
if(size(m,2) == 1) %��mת����������
    m = m';
end
m = t * m;%��������mת���ɸ�����ͬ�ľ�������ͬk����
value = ones(size(m));
t = value;
for ii = 0 : L
    tt = ii * t;
    k_tt = k - tt;
    m_tt = m - tt;
    [i,jj,k_tt_0] = find(abs(k_tt) < eps);
    [mm,nn]=size(k_tt);
    k_tt_0 = sparse(i,jj,k_tt_0,mm,nn);
    k_tt_0 = full(k_tt_0);
    k_tt_0 = m_tt .* k_tt_0;
    
    k_tt = k_tt + k_tt_0;
    
    value = value .* m_tt ./ k_tt;
end

function rt2 = f(rt)

function value = bmk(j,k,m,L)
% ������kת���ɸ�����ͬ�ľ���������m��ͬ
t = ones(1,length(m));
if(size(k,1) == 1) %��kת����������
    k = k';
end
k = k * t;
%��������mת���ɸ�����ͬ�ľ�������ͬk����
t = ones(size(k,1),1);
if(size(m,2) == 1) %��mת����������
    m = m';
end
m = t * m;%��������mת���ɸ�����ͬ�ľ�������ͬk����
value = ones(size(m));
t = value;
for ii = (2.^j - L) : (2.^j)
    tt = ii * t;
    k_tt = k - tt;
    m_tt = m - tt;
    [i,jj,k_tt_0] = find(abs(k_tt) < eps);
    [mm,nn]=size(k_tt);
    k_tt_0 = sparse(i,jj,k_tt_0,mm,nn);
    k_tt_0 = full(k_tt_0);
    k_tt_0 = m_tt .* k_tt_0;
    
    k_tt = k_tt + k_tt_0;
    
    value = value .* m_tt ./ k_tt;
end
%========================================================================
function [d0theta,d1theta,d2theta]=dxfai0(j,k,r,xmin,xmax,x) %theta��һ�׺Ͷ��׵���ֵ
%����k�����������
%����x�����������
%d1theta,d2theta���Ǿ���������kȷ����������x��ά��ȷ��
global delta
global mm

delta = (xmax-xmin)/2.^j;
mm = r*r*delta*delta;

xk = xmin + delta * k;  % xk������

%������xת���ɸ�����ͬ�ľ���������k
t = ones(1,length(k));
if(size(x,1) == 1) %��xת����������
    x = x';
end
x = x * t; %��������xת���ɸ�����ͬ�ľ���������k

%��������xkת���ɸ�����ͬ�ľ�������ͬx����
t = ones(size(x,1),1);
if(size(xk,2) == 1) %��xkת����������
    xk = xk';
end
xk = t * xk;%��������xkת���ɸ�����ͬ�ľ�������ͬx����

deltaxk = x - xk;  %deltaxk�Ǿ���

[i,j,a] = find(abs(deltaxk) < eps);
[m,n]=size(deltaxk);
a = sparse(i,j,a,m,n);
a = full(a);
d0theta = a;
d1theta = 0 * a;
d2theta = -(3+pi*pi*r*r)/(3*mm) * a;

c0 = spfun(@sad0x,deltaxk);
c1 = spfun(@sad1x,deltaxk);%���պ���sada����deltaxk�з�0Ԫ��
c2 = spfun(@sad2x,deltaxk);%���պ���sada����deltaxk�з�0Ԫ��
c0 = full(c0);
c1 = full(c1);
c2 = full(c2);
d0theta = d0theta + c0;
d1theta = d1theta + c1;
d2theta = d2theta + c2;
%---------------------------------------
function d0theta = sad0x(deltaxk)
global delta
global mm

s = sin(pi * deltaxk / delta);
e = exp(-(deltaxk .^ 2) / (2 * mm));

se = s .* e;

d0theta = se ./ (pi * deltaxk ./ delta);
%---------------------------------------
function d1theta = sad1x(deltaxk)
global delta
global mm

s=sin(pi*deltaxk/delta);
c=cos(pi*deltaxk/delta);
e=exp(-(deltaxk.^2)/(2*mm));

se=s.*e;
ce=c.*e;

d1theta=ce./deltaxk-se.*delta./(pi*deltaxk.^2)-se.*delta/(pi*mm);
%---------------------------------------
function d2theta = sad2x(deltaxk)
global delta
global mm

s=sin(pi*deltaxk/delta);
c=cos(pi*deltaxk/delta);
e=exp(-(deltaxk.^2)/(2*mm));

se=s.*e;
ce=c.*e;

d2theta=-pi*se./(delta*deltaxk)-2*ce./(deltaxk.^2)-...
    2*ce/mm+2*se*delta./(pi*deltaxk.^3)+se*delta./(pi*mm*deltaxk)+...
    deltaxk.*se*delta/(pi*mm*mm);

%========================================================================
%����С��fai(x-xn)��ֵ��������
function [value,dvalue,ddvalue] = interval012d(j,r,xmin,xmax,x)
%����x��xnum�����������,xȷ������,xnumȷ������
L = 1;%the number of extern wavelet collocation point
N = 2 * j;

xnum = 0 : 2 ^ j;
[value,dvalue,ddvalue] = dxfai0(j,xnum,r,xmin,xmax,x);
deltax = (xmax-xmin) / 2 ^ j;
s = 0; ds = 0; dds = 0;

xnum = 0 : L;
n =  -N + 1 : -1;
ankvalue = amk(j,xnum,n,L);
[Ii,Ii1d,Ii2d] = dxfai0(j,n,r,xmin,xmax,x);
s = Ii * ankvalue';
ds = Ii1d * ankvalue';
dds = Ii2d * ankvalue';

value(:,xnum+1) = value(:,xnum+1) + s;
dvalue(:,xnum+1) = dvalue(:,xnum+1) + ds;
ddvalue(:,xnum+1) = ddvalue(:,xnum+1) + dds;

xnum = 2^j - L : 2 ^ j;
n = (2 ^ j + 1) : (2 ^ j + N - 1);
bnkvalue = bmk(j,xnum,n,L);
[Ii,Ii1d,Ii2d] = dxfai0(j,n,r,xmin,xmax,x);
s = Ii * bnkvalue';
ds = Ii1d * bnkvalue';
dds = Ii2d * bnkvalue';

value(:,xnum+1) = value(:,xnum+1) + s;
dvalue(:,xnum+1) = dvalue(:,xnum+1) + ds;
ddvalue(:,xnum+1) = ddvalue(:,xnum+1) + dds;
