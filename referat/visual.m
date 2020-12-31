clear;
n = 50; % Number of space steps
L = 1; % Length of the wire
T =20; % Final time
maxk = 300000; % Number of time steps
dt = T/maxk;
dx = (L)/n;
b = dt/(dx*dx); % Stability parameter (b=<1)
for i = 1:n+1
    if(mod(i,5)==0)
        f(i) =50;
    else
        f(i)=0;
    end
end
for i = 1:n+1
    if(mod(i,6)==0)
        f(i) =50;
    else
        f(i)=0;
    end
end
for i = 1:n+1
    if(mod(i,7)==0)
        f(i) =50;
    else
        f(i)=0;
    end
end
for i = 1:n+1
x(i) =(i-1)*dx;
u(i,1) =24;
end
for k=1:maxk+1
time(k) = (k-1)*dt;
end

for k=1:maxk % Time Loop
for i=2:n % Space Loop
u(i,k+1) =u(i,k) + 0.001*b*(u(i-1,k)+u(i+1,k)-2.*u(i,k))+0.5*dt*f(i);
end
u(1,k+1)=4/3*u(2,k+1)-1/3*u(3,k+1);
u(n+1,k+1) =4/3*u(n,k+1)-1/3*u(n-1,k+1);
end
% Graphical representation of the temperature at different selected times
figure(1)
grid on
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,1000),'-',x,u(:,2000),'-')
title('Temperature analysis in time moments t=0.001, t=0.1, t=0.3, t=1, t=2')
xlabel('X')
ylabel('T')
figure(2)
grid on
mesh(x,time,u')
title('Visual representation of the temperature distribution on the plate')
xlabel('X')
ylabel('Temp')
