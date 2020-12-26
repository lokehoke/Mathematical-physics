clc; clear; close all;
c=10;
syms x;
syms t;
syms tau;
pi=3.14;
l=1; % Длина участка струны
ro=0.03; % Плотность
a=6; %амплитуда
omega=0.8; %частота
save=0;
 % U_tt-c^2*U_xx=f(t,x)
u_0x=sin(3.14*x);
Ut_0x=0;
u_t0=0;
u_tl=0;
st=1;
N=15;
h=0.008; %шаг по координате
tu=0.000625; %шаг по времени
X0 = 0; Xl = l; % Границы по X
T0 = 0; Tl = 5; % Границы по T
M=round((Xl-X0)/h); %число разбиений по координате
P=round((Tl-T0)/tu); %число разбиений по времени
U1=0;
U2=0;
st=5;
for n=1:N
    k(n)=pi*n/l;
    omeg(n)=pi*n*c/l;
    phi(n)=2/l*int(u_0x*sin(k(n)*x),x,0,l);
    xi(n)=2/l*int(Ut_0x*sin(k(n)*x),x,0,l);
    if (rem(n,2)==0)
       T(n)=0;
    else
       T(n)=((4*a)/(omeg(n)*pi*l*ro*n))*int(sin(omega*tau)*sin(omeg(n)*(t-tau)),tau,0,t);
    end
    U1=U1+(phi(n)*cos(omeg(n)*t)+xi(n)/phi(n)*sin(omeg(n)*t))*sin(k(n)*x);
    U2=U2+T(n)*sin(k(n)*x);
end
U=U1+U2;
 
%Q=zeros(M+1,P+1);
 
for p=1:300
   Q(p)=subs(U,t,T0+tu*(p-1));
end
 
figure;
if(save==1)
hold on;
end
for i=1:1:300
    A=Q(i);
ezplot(A,X0,Xl);
axis([0 l -2 2]); 
grid on;
xlabel('x'); ylabel('U(t,x)');
title(sprintf('Колебания струны при t = %fс',T0+tu*(i-1)));
pause(0.01)
end
