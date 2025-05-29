function [t,theta]=laneemdenode45(n,epsilon,puntofinal)
%Resuelve la ecuación de Lane-Emden para un índice politrópico n mediante
%ode45
%Entradas:
%   n: Índice politrópico (puede ser cualquier valor)
%   epsilon: Valor inicial pequeño para evitar la singularidad en t=0
%   puntofinal: Valor máximo de t hasta donde integrar
%Salidas:
%   theta: Solución de la ecuación de Lane-Emden (theta(t))
%   t: Puntos de la coordenada donde se calcula la solución

%Condiciones iniciales en t=epsilon
theta0=1; %Aproximación
dtheta0=0; %Derivada de la aproximación
%Definir el sistema de ecuaciones para ode45
funcion=@(t,y)[y(2);-y(1)^n-(2/t)*y(2)];
%Resolver desde epsilon hasta puntofinal
[t,y]=ode45(funcion,[epsilon,puntofinal],[theta0;dtheta0]);
%Extraer theta de la solución
theta=y(:,1);
%Si se alcanza theta=0 (es decir la solucion) antes del puntofinal
%la solución la mostraraá haciendo lo siguiente, que es cuando cruza 
zero_crossing=find(theta <= 0, 1); 
if ~isempty(zero_crossing)           
    t=t(zero_crossing(end));
    theta=theta(zero_crossing(end));
end
