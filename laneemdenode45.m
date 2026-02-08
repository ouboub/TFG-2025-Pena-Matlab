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
% Esta función resuelve la ecuación de Lane-Emden para un índice
% politrópico n mediante ode45.
% Entradas:
%   n: Índice politrópico.
%   epsilon: Valor inicial pequeño para evitar la singularidad en t=0.
%   puntofinal: Valor máximo de t hasta donde se puede integrar.
% Salidas:
%   theta: Solución de la ecuación de Lane-Emden (theta(t)).
%   t: Primer cero de la ecuación.

% Condiciones iniciales en t = epsilon.
theta0=1;  % Aproximación.
dtheta0=0; % Derivada de la aproximación.
% Definimos el sistema de ecuaciones para ode45.
funcion=@(t,y)[y(2);-y(1)^n-(2/t)*y(2)];
% Resolvemos desde epsilon hasta puntofinal.
[t,y]=ode45(funcion,[epsilon,puntofinal],[theta0;dtheta0]);
% Extraemos theta de la solución.
theta=y(:,1);
% Si se alcanza theta=0 (es decir, la solución) antes del punto final del
% intervalo, el primer cero se mostrará haciendo lo siguiente. 
% Aquí me he ayudado de ChatGPT, porque no sabía como podía programar que
% me diera el primer cero cuando cruza.
zero_crossing=find(theta <= 0, 1); 
if ~isempty(zero_crossing)           
    t=t(zero_crossing(end));
    theta=theta(zero_crossing(end));
end
% Para poder dibujar la gráfica necesitariamos comentar estas últimas cinco
% líneas, ya que sino solo dibujaría un punto, el de la solución.