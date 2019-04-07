function [xlagrange, ylagrange, zlagrange] = GetLagrangePoints(N)
global M_earth M_moon

CI = [0.5,0,0;
    1.5,0,0;
    -0.5,0,0;
    0.5,0.5,0;
    0.5,-0.5,0];    %condiciones iniciales, se eligen en función del punto de Lagrange buscado


[xyz]=fsolve(@(U) LagrangeCoordinates(U),CI(N,:));

xlagrange = xyz(1);
ylagrange = xyz(2);
zlagrange = xyz(3);    


function [U_sol] = LagrangeCoordinates(U)

M_earth = 5.972e24;  %kg
M_moon = 7.349e22;   %kg
mu = M_moon / (M_moon + M_earth);

U4=0;
U5=0;
X = U(1) + 2*U5- ( 1-mu)*(U(1)+mu )/(sqrt( (U(1)+mu)^2 + U(2)^2 + U(3)^2))^3 - mu*(U(1)-1+mu)/(sqrt( (U(1)-1+mu)^2 + U(2)^2 + U(3)^2))^3;
Y = U(2) + 2*U4 - (1-mu)*U(2)/(sqrt( (U(1)+mu)^2 + U(2)^2 + U(3)^2))^3 - mu*U(2)/(sqrt( (U(1)-1+mu)^2 + U(2)^2 + U(3)^2))^3;
Z = - (1-mu)*U(3)/(sqrt( (U(1)+mu)^2 + U(2)^2 + U(3)^2))^3 - mu*U(3)/(sqrt( (U(1)-1+mu)^2 + U(2)^2 + U(3)^2))^3;

U_sol = [X Y Z];
end


end