%   OSCILADOR ARMÓNICO CON EULER EXPLÍCITO

clc
clear all
close all
format shortg

x0=0;
v0=1;
t_inicial=0;
t_final = 100;    %1 año en segundos
N=1e4;

U0 = [x0;v0];
A=[0,1;-1,0];


U_oscilador = ProblemaCauchy(U0,t_inicial,t_final,N,A,@EulerInverso);
plot(U_oscilador(1,:),U_oscilador(2,:))
title('Oscilador armónico Euler inv.')
xlabel('x')
ylabel('dxdt')

function [U] = ProblemaCauchy(U0,t1,t2,N,F,Integrador)    %Wrapper problema Cauchy
    Delta_t = (t2-t1)/N;    %Paso temporal constante
    U(:,1) = U0;    %condición inicial
    for i=2:N
        U(:,i) = Integrador(U(:,i-1),Delta_t,F);
    end
end

function [U] = EulerInverso(U0,Dt,F)

    U = (eye(size(F)) - Dt.*F)\U0;
end
 
function [F] = OpDifOscilador(U)
    F = [U(2);-U(1)];
end
    
