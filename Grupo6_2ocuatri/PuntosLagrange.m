function [U_Lagrange,Autoval,Error] = PuntosLagrange(U0,t1,t2,N,Integrador,q)
global M_earth M_moon

%   MATRIZ JACOBIANA Y ESTABILIDAD
epsilon = zeros(length(U0),1);

for k=1:length(U0)
    epsilon(:)=0;
    epsilon(k) = 1e-3;
    Jacobiano(:,k) = (Orbitas(U0' + epsilon(:)) - Orbitas(U0' - epsilon(:)))  ./ norm(2*epsilon(:));
end

[autovec, autoval] = eig(Jacobiano);
Autoval = diag(autoval);

%%   PROBLEMA CAUCHY
for i=1:2
    Time_domain = linspace(t1,t2,i*N);
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    if i==1
        U(1,:) = U0;    %condición inicial
        [t,U_sol1] = Integrador( @(t,U)Orbitas(U), Time_domain, U0,options); 
    else
        U0 = U0+0.001;
        U(1,:) = U0;
        [t,U_sol2] = Integrador( @(t,U)Orbitas(U), Time_domain, U0,options); 
    end
end
    U_Lagrange = U_sol1;
    
    
%%   EXTRAPOLACIÓN DE RICHARDSON

dt = (t2-t1)/N;

K_Richardson = zeros(N,1);
for i=1:N
K_Richardson(i) = norm(U_sol2(2*i-1,:) - U_sol1(i,:))./(-(dt/2)^q+dt^q);
end

Error = dt^q.*K_Richardson;


%%   FUNCIÓN DEL SISTEMA EN EJES NO INERCIALES
function F = Orbitas(U)
    mu = M_moon / (M_moon + M_earth);
    d = sqrt( (U(1)+mu)^2 + U(2)^2 + U(3)^2);
    r = sqrt( (U(1)-1+mu)^2 + U(2)^2 + U(3)^2);
    
    F=zeros(1,6);
    F(1:3) = U(4:6);
    F(4) = U(1) + 2*U(5)- ( 1-mu)*(U(1)+mu )/d^3 - mu*(U(1)-1+mu)/r^3;
    F(5) = U(2) - 2*U(4) - (1-mu)*U(2)/d^3 - mu*U(2)/r^3;
    F(6) = - (1-mu)*U(3)/d^3 - mu*U(3)/r^3;
    
    F=F';
end

end