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

function [alfa,y,alfas,a,x]=laneemdenNEWTONKANTOROVICH(n,N,maxiter)
    %   Resuelve la ecuación de Lane-Emden usando el método espectral de Chebyshev
    %   y el de Newton-Kantorovich
    %   ENTRADAS:
    %   n: exponente en la ecuación Lane-Emden 
    %   N: número de puntos de colocación
    %   maxiter: máximo número de iteraciones de Newton 
    %   SALIDAS:
    %   alfa: autovalor (solución de la ecuación) 
    %   y: solución evaluada en x 
    %   alfas: distintas soluciones para cada n
    %   x: puntos de la malla
    %   a: coeficientes de chebyshev

    % 1.Puntos de Chebyshev-Lobatto en [0,1]
    j=(0:N-1)';
    tj=pi*(j/(N-1));  % Puntos en [-1,1]
    x=(1+cos(tj))/2;  % Transformacion a [0,1]    
    % 2.Matrices de diferenciación de Chebyshev
    [D0,D1,D2]=chebdiff(N);
    % 3.Condiciones iniciales
    y0=cos((pi/2).*x); 
    alfa=3;
    % 4.Encontrar coeficientes iniciales de Chebyshev  
    a=D0\y0;
    % 5.Iteraciones de Newton-Kantorovich
    alfas=zeros(maxiter,1);
    alfas(1)=alfa;
    for iter=1:maxiter
        % Evaluaamos y y sus derivadas
        y=D0*a;
        yx=D1*a;
        yxx=D2*a;
        % Iniciamos el sistema lineal para la corrección (Jacobiano)
        J=zeros(N+1,N+1);
        r=zeros(N+1,1);
        % Ecuación en puntos interiores (N-2 ecuaciones)
        J(1:N-2,1:N)=D2(2:N-1,:)+(2./x(2:N-1)).*D1(2:N-1,:)+alfa^2*n*y(2:N-1).^(n-1).*D0(2:N-1,:);
        J(1:N-2,N+1)=2*alfa*y(2:N-1).^n; 
        r(1:N-2)=yxx(2:N-1)+(2./x(2:N-1)).*yx(2:N-1)+alfa^2.*y(2:N-1).^n;
        % Condiciones de frontera
        % y(0)=1
        J(N-1,1:N)=D0(1,:);
        r(N+1)=y(N)-1;
        % yx(0)=0 
        J(N,1:N)=D1(end,:);
        % y(1)=0
        J(N+1,1:N)=D0(end,:);
        % Resolver el sistema lineal
        delta=-J\r;
        % Actualizar solución
        a=a+delta(1:N);
        alfa=alfa+delta(end);
        alfas(iter)=sum(alfa);
    end
    % 6.Evaluar la solucion final
    y=D0*a;
end
function [D0,D1,D2]=chebdiff(N)
    % Construye matrices de diferenciación de Chebyshev para [0,1]
    j=0:N-1;
    theta=pi*j/(N-1);
    t=theta(:);  % columna
    J=j(:)';     % fila
    T=cos(J.*t);         % Matriz D0(NxN)
    S=sin(t); 
    C=cos(t);
    D0=T;
    % Construimos D1 y D2 usando fórmulas
    D1=2*(J.*sin(J.*t))./S;
    D2=4*(-(J.^2.*cos(J.*t))./(S.^2)+(J.*sin(J.*t).*C)./(S.^3));
    % Corregir extremos
    for k=1:N
        n=k-1;
        D0(1,k)=1;
        D0(N,k)=(-1)^n;
        D1(1,k)=2*n^2;
        D1(N,k)=2*(-1)^n*(n^2);
        D2(1,k)=4/3*(n^2)*(n^2-1);
        D2(N,k)=4/3*(-1)^n*(n^2)*(n^2-1);
    end
end

