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
function [alfa,y,alfas]=laneemdenSLM(n,N,maxiter)
% Función que resuelve la ecuación de Lane-Emden para un N mayor que 2 y un
% n comprendido entre 0 y 5. 
% ENTRADAS
    %   n = índice politrópico de la ecuación.
    %   N = número de puntos de colocación.
    %   maxiter = máximo número de iteraciones.
% SALIDAS
    %   alfa = parámetro (primer cero de la ecuación).
    %   y = solución aproximada.
    %   alfas = parámetros para cada n.

    % Puntos iniciales de Gauss-Lobato.
    z=cos(pi*(0:N)'/N);
    % Tansformación del dominio.
    t=(z+1)/2;
    % Matriz de Chebyshev, más abajo definida.
    D=chebdiff(N);
    D1=2*D;
    D2=4*D^2;
    % Evaluación inicial de alfa.
    alfa0=sqrt(6*(4/3)^n);
    % Valor inicial de y.
    y0=1-t.^2;
    % Almacenamiento de iteraciones.
    ytotal=y0;
    alfatotal=alfa0;
    alfas=zeros(maxiter,1);
    alfas(1)=alfa0;
    % Iteraciones del SLM.
    for iter=1:maxiter
        % Sumamos las componentes mediante el comando sum.
        suma_y=ytotal;
        suma_alfa=alfatotal; 
        % Calculamos las matrices a1 y a2.
        a1=n*suma_alfa^2*suma_y.^(n-1);
        a2=2*suma_alfa*suma_y.^n;
        % Matrices de diferenciación para aproximar las derivadas.
        dy=D1*suma_y;
        ddy=D2*suma_y;
        % Cálculo del residuo.
        r=-(ddy+(2./t).*dy+suma_alfa^2*suma_y.^n);
        % Construimos la matriz A y B.
        A=D2+diag(2./t)*D1+diag(a1);
        B=a2;
        % Eliminación de la primera y última fila/columna.
        Amatriz=A(2:N,2:N);
        Bmatriz=B(2:N);
        Rmatriz=r(2:N);
        % Definimos las matrices aumentadas.
        Asistema=[Amatriz,Bmatriz;D(N+1,2:N),0];
        Rsistema=[Rmatriz;0];
        % Resolvemos el sistema aumentado de todas las matrices juntas.
        solucion=Asistema\Rsistema;
        % Extraemos la solución una a una del sistema aumentado.
        Y=zeros(N+1,1);
        Y(2:N)=solucion(1:end-1);
        alfa_i=solucion(end);
        % Actualización de los valores.
        ytotal=ytotal+Y;
        alfatotal=alfatotal+alfa_i;
        alfas(iter)=sum(alfatotal);
    end
    % Soluciones finales.
    alfa=alfatotal;
    y=ytotal;
end
function D=chebdiff(N)
% CHEBDIFF. Matriz de diferenciación de Chebyshev.
% D=CHEBDIFF(N) calcula la matriz de diferenciación de Chebyshev de tamaño (N+1)x(N+1).

% Puntos de Chebyshev-Gauss-Lobatto.
z=cos(pi*(0:N)'/N);
% Coeficientes de la matriz de diferenciación de Chebyshev.
c=[2;ones(N-1,1);2];
c=c.*(-1).^(0:N)';
% Matriz de diferenciación, elementos D(j,K).
X=repmat(z,1,N+1);
dX=X-X';
D=(c*(1./c)')./dX;
D(eye(N+1)~=0)=0; % Se anulan los elementos diagonales.
% Nuevos elementos diagonales.
for k=1:N-1
    D(k+1,k+1)=-z(k+1)/(2*(1-z(k+1)^2));
end
% Condiciones para 00 y NN.
D(1,1)=(2*N^2+1)/6;
D(N+1,N+1)=-D(1,1);
end
  