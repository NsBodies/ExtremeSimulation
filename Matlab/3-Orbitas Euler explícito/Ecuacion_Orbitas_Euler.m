%ECUACION DE LAS ORBITAS (EULER EXPLICITO)

clc
clear all
close all
format shortg

x0=[0;6900];
v0=[7.6005;0];  %Velocidad para órbita circular

t_inicial=0;
t_final = 365*24*3600;    %1 día en segundos
N=1e6;

U0 = [x0;v0];

U_orbitas = ProblemaCauchy(U0,t_inicial,t_final,N,@OpDifOrbitas,@EulerExplicito);
plot(U_orbitas(1,:),U_orbitas(2,:))
title('Órbita Euler expl. (1 día)')
xlabel('x (km)')
ylabel('y (km)')

function [U] = ProblemaCauchy(U0,t1,t2,N,F,Integrador)    %Wrapper problema Cauchy
    Delta_t = (t2-t1)/N;    %Paso temporal constante
    U(:,1) = U0;    %condición inicial
    for i=2:N
        U(:,i) = Integrador(U(:,i-1),Delta_t,F);
    end
end

function [U] = EulerExplicito(U0,Dt,F)
F_paso = F(U0);
    U = U0 + Dt.*F_paso;
end
 
function [F] = OpDifOrbitas(U)
    GM = 3.986e5;  
    R = norm(U(1:2));
    F = [U(3);U(4);-GM./R^2.*(U(1)/R);-GM./R^2.*(U(2)/R)];

end
    

