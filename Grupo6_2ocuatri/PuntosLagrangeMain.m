%   PUNTOS LAGRANGE GUAPOS

clc
clear all
close all
format shortg


global M_earth M_moon
    M_earth = 5.972e24;  %kg
    M_moon = 7.349e22;   %kg

P = 3; %INDICAR AQUÍ EL PUNTO DE LAGRANGE BUSCADO

[r0(1), r0(2), r0(3)] = GetLagrangePoints(P);  
v0 = [0.01, 0, 0];    %Aquí se incluye la perturbación
U0 = [r0,v0];   %Vector de condiciones iniciales

t_inicial = 0;
t_final = 100;
N=1e5;  %Número de pasos temporales
q=4;    %Orden del esquema numérico a utilizar

[U_orbitas, Autoval,Error] = PuntosLagrange(U0,t_inicial,t_final,N,@ode45,q); %RK Orden 4
[U_orbitas_ABM, Autoval,Error] = PuntosLagrange(U0,t_inicial,t_final,N,@ode113,q);  %ABM






%   PLOTEADO
figure
plot(U_orbitas(:,1),U_orbitas(:,2))
title('Órbita mediante Runge-Kutta de orden 4');
xlabel('$$\hat{x}$$','Interpreter','Latex');ylabel('$$\hat{y}$$','Interpreter','Latex')
figure
plot(U_orbitas_ABM(:,1),U_orbitas_ABM(:,2),'r')
title('Órbita mediante Adam-Bashford-Moulton')
xlabel('$$\hat{x}$$','Interpreter','Latex');ylabel('$$\hat{y}$$','Interpreter','Latex')

Diff = max(abs(U_orbitas_ABM-U_orbitas)); %Máxima difrencia entre los métodos Runge-Kutta y ABM

figure
plot(linspace(t_inicial,t_final,N),abs(U_orbitas_ABM-U_orbitas))
title('Difrencia entre los métodos Runge-Kutta y ABM para cada paso temporal');
xlabel('t'); ylabel('Diferencia');
legend('X','Y','Z','V_x','V_y','V_z')

% figure
% plot(linspace(t_inicial,t_final,N),Error)


