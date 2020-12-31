clc; clear; close all;
c=10;
l=1; % Длина участка струны
ro=0.03; % Плотность
a=6; %амплитуда
omega=0.8; %частота
f=@(t,x)(a*sin(omega*t))/(l*ro); % U_tt-c^2*U_xx=f(t,x)
u_0x=@(x)sin(3.14*x);
Ut_0x=@(x)0;
u_t0=@(x)0;
u_tl=@(x)0;
st=1;
 
save=0; % Способ вывода: 0 - без сохранения, 1 - с сохранением
% Расчет шагов сетки
X0 = 0; Xl = l; % Границы по X
T0 = 0; Tl = 5; % Границы по T
N=1000; % Количество разбиений по X
h=(Xl-X0)/N; tau=(Tl-T0)/N;
while(c*tau/h>=1) % Проверка сходимости метода
tau=tau/2;
h=h*2;
end
M=(Xl-X0)/h;
P=(Tl-T0)/tau;
% Создание сетки
U=zeros(M+1,P+1);
 
for m=1:M+1 %U(m,0)(m = 0, 1, ..., M) и U(m,1)(m = 1, 2, ..., M ? 1)
    U(m,1)=u_0x(X0+(m-1)*h);
    U(m,2)=tau*Ut_0x(X0+(m-1)*h)+U(m,1);
end
 
%U(0,1) и U(M,1)
 U(1,2)=u_t0(T0+tau);
 U(M+1,2)=u_tl(T0+tau);
 
 
for p=2:P
    for m=2:M
    U(m,p+1)=(tau^2)*(f(T0+tau*(m-1),X0+h*(p-1))+((c^2)/(h^2))*(U(m+1,p)-2*U(m,p)+U(m-1,p)))+2*U(m,p)-U(m,p-1);
    end
end
 
if(max(max(abs(U)))>10^6)
error('Расхождение схемы')
end
 
X = X0:h:Xl;
figure;
if(save==1)
hold on;
end
for i=1:st:P+1
plot(X,U(:,i),'b-','linewidth',2);
axis([0 l -1.4 1.4]); 
grid on;
xlabel('x'); ylabel('U(t,x)');
title(sprintf('Колебания струны при t = %fс',T0+tau*(i-1)));
pause(0.001)
end
