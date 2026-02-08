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
function [t,theta]=laneemdenRK4(n,epsilon,puntofinal,N)
    % Método de Runge-Kutta de orden 4 para resolver la ecuación de
    % Lane-Emden.
    % ENTRADAS:
    %   n: Índice politrópico.
    %   epsilon: Valor inicial pequeño.
    %   puntofinal: Punto final de integración.
    %   N: Paso fijo para RK4.
    % SALIDAS:
    %   t: Primer cero de la ecuación.
    %   theta: Solución de la ecuación evaluada en t.

    % Inicialización
    t=epsilon:N:puntofinal; % Vector de puntos de la malla. 
    M=length(t);  
    y=zeros(M,2); % Matriz de la solución.
    theta0=1;
    dtheta0=0;
    y(1,:)=[theta0,dtheta0]; % Condiciones iniciales.
    funcion=@(t,y)[y(2);-(2/t)*y(2)-y(1)^n]; % Ecuación de Lane-Emden.
    % Iteración RK4.
    for i=1:M-1
        ti=t(i);
        yi=y(i,:)';
        % Etapas RK4.
        k1=funcion(ti,yi);
        k2=funcion(ti+N/2,yi+(N/2)*k1);
        k3=funcion(ti+N/2,yi+(N/2)*k2);
        k4=funcion(ti+N,yi+N*k3);
        % Actualización.
        y(i+1,:)=yi+(N/6)*(k1+2*k2+2*k3+k4);
    end
% Sacamos la solución.
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