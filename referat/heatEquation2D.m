clear;

plane_num = 52;
tube_num = 8;

p = 100;
T_simulation = 10; % time of simulation
T = 5000; % num of dt
NP = 5; % num of dx in plane
NT = 5; % num of dx in tube

lp = 0.458; % length of plane
lt = 0.13; % length of tube

dt = T_simulation / T;
dx_p = lp / NP;
dx_t = lt / NT;
d_tube = 0.006;
dt / dx_p / dx_p
dt / dx_t / dx_t

for (j = 1:plane_num)
    for (n = 1:NP)
        U(j, n, 1) = 24;
    end
end


for (i = 1:tube_num)
    for (n = 1:NT)
        V(i, n, 1) = 24;
    end
end

f_tube(1) = 0;
f_plate(1) = 0;

width_plane = 0.00035;
diff_plane = 0.0002; % distance between plates


% find field where plane touch tube
j = 1;
for (i = 0:dx_t:lt)
    target_tube(j) = 0;

    current_plane = mod(i, width_plane + diff_plane);
    if (current_plane < width_plane) 
        target_tube(j) = 1;
    end

    j = j + 1;
end

% find field where tube touch plane 
x_center = [0.014, 0.042, 0.088, 0.114];
j = 1;
for (i = 0:dx_p:lp)
    target_plane(j) = 0;    
    for (k = 1:length(x_center))
        if (i > x_center(k) - d_tube / 2 && i < x_center(k) + d_tube / 2)
            target_plane(j) = 1;
        end
    end

    j = j + 1;
end

k_diff = 2;
F_cooler = 1/100000;
rpm = 1400;

for (t = 1:T-1) % Time
    t
    V_mean = mean(mean(V, 1), 2);
    U_mean = mean(mean(U, 1), 2);

    %%% calculate tube
    for (i = 1:tube_num)
        for (n = 2:NT-1)
            if (t ~= 1) 
                f_tube(t) = -k_diff*(V_mean(t) - U_mean(t-1))*target_tube(n);
            end
            V(i, n, t+1) = V(i, n, t) + 2*dt/dx_p^2*(V(i, n+1, t) - 2*V(i, n, t) + V(i, n-1, t)) + dt*f_tube(t);
        end
        V(i, 1, t+1) = 4/3*V(i, 2, t+1) - 1/3*V(i, 3, t+1);
        V(i, NT, t+1) = 4/3*V(i, NT-1, t+1) - 1/3*V(i, NT-2,t+1) + 2/3*dx_t*p/4;
    end
   
    
    V_mean = mean(mean(V, 1), 2);
    U_mean = mean(mean(U, 1), 2);
    
    %%% calculate plane
    for (j = 1:plane_num)
        for (n = 2:NP-1)
            if (t ~= 1) 
                f_plate(t) = k_diff*(V_mean(t) - U_mean(t-1))*target_plane(n) + F_cooler*rpm;
            end
            U(j, n, t+1) = U(j, n, t) + 2*dt/dx_p^2*(U(j, n+1, t) - 2*U(j, n, t) + U(j, n-1, t)) + dt*f_plate(t);
        end
        U(j, 1, t+1) = 4/3*U(j, 2, t+1) - 1/3*U(j, 3, t+1);
        U(j, NP, t+1) = 4/3*U(j, NP-1, t+1) - 1/3*U(j, NP-2,t+1);
    end

    
end

U_mean = mean(mean(U, 1), 2);
V_mean = mean(mean(V, 1), 2);

figure(1);
plot(0:dt:T_simulation-dt, U_mean);
title('mean temperature of planes')
xlabel('t, c');
ylabel('T, C');
colormap winter;

figure(2);
plot(0:dt:T_simulation-dt, V_mean);
title('mean temperature of tubes')
xlabel('t, c');
ylabel('T, C');
colormap winter;

figure(3);
[X, Y] = meshgrid(0:dt:T_simulation-dt, 0:dx_p:lp-dx_p);
mesh(X, Y, squeeze(mean(U, 1)));
title('mean temperature of planes in each point')
xlabel('t, c');
ylabel('x, m');
zlabel('T, C');
colormap winter;


figure(4);
[X, Y] = meshgrid(0:dt:T_simulation-dt, 0:dx_t:lt-dx_t);
mesh(X, Y, squeeze(mean(V, 1)));
title('mean temperature of planes in each point')
xlabel('t, c');
ylabel('x, m');
zlabel('T, C');
colormap winter;
