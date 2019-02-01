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

U_oscilador = ProblemaCauchy(U0,t_inicial,t_final,N,@OpDifOscilador,@EulerExplicito);
plot(U_oscilador(1,:),U_oscilador(2,:))
title('Oscilador armónico Euler expl.')
xlabel('x')
ylabel('dxdt')

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
 
function [F] = OpDifOscilador(U)
    F = [U(2);-U(1)];
end
    
