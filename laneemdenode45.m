% Copyright (C) 2026 Diego Fernandez de la Peña

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
%la solución la mostraraá haciendo lo siguiente, que es cuando cruza.
%aquí me he ayudado de chatgpt porque no sabia como podia programar que la
%me diera la solucion cuando cruza
zero_crossing=find(theta <= 0, 1); 
if ~isempty(zero_crossing)           
    t=t(zero_crossing(end));
    theta=theta(zero_crossing(end));
end
%para poder dibujar la grafica necesitariamos comentar estas últimas cinco
%líneas ya que si no solo dibujaria un punto, el de la solución