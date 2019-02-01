%ECUACION DE LAS ORBITAS (ODE45)

clc
clear all
close all
format shortg

x0=[0;6900];    %Posicion inicial (km)
v0=[7.6005;0];  %Velocidad para órbita circular (km/s)

t_inicial = 0;
t_final = 365*24*3600;    %1 día en segundos
N = 1e5;

U0 = [x0;v0];



U_orbitas = ProblemaCauchyRunge(U0,t_inicial,t_final,N,@OpDifOrbitas,@ode45);
plot(U_orbitas(1,:),U_orbitas(2,:))
title('Órbita RK4 (1 año)')
xlabel('x (km)')
ylabel('y (km)')
grid on

function [U] = ProblemaCauchyRunge(U0,t1,t2,N,F,Integrador)    %Wrapper problema Cauchy
    U(:,1) = U0;    %condición inicial
    Time_domain = linspace(t1,t2,N);
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    [t,U] = Integrador( @(t,U)F(U), Time_domain, U0,options);
     U = U';
end

 
function F = OpDifOrbitas(U)
    GM = 3.986e5;  
    R = norm(U(1:2));
    F = [U(3);U(4);-GM/R^2*(U(1)/R);-GM/R^2*(U(2)/R)];
end
    
