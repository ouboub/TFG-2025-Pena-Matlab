function [alfa,y,alfas]=laneemdenSLM(n,N,maxiter)
% Funcion que resuelve la ecuacion de Lane-Emden para un N mayor que 2 y un
% n comprendido entre 0 y 5 
% ENTRADAS
    %   n=índice politrópico de la ecuación
    %   N=número de puntos de colocación (tamaño)
    %   maxiter=máximo número de iteraciones
% SALIDA
    %   alfa=valor propio (primera solucion de y)
    %   y=solución aproximada
    %   alfas=distintas soluciones para cada n
    %Puntos iniciales de Gauss-Lobato
    z=cos(pi*(0:N)'/N);
    %Tansformacion del dominio
    t=(z+1)/2;
    %Matriz derivada de chebyshev más abajo definida
    D=chebdiff(N);
    D1=2*D;
    D2=4*D^2;
    %Evaluacion inicial de alfa
    alfa0=sqrt(6*(4/3)^n);
    %Valor inicial de y
    y0=1-t.^2;
    %Almacenamiento de iteraciones
    ytotal=y0;
    alfatotal=alfa0;
    alfas=zeros(maxiter,1);
    alfas(1)=alfa0;
    %Iteraciones del SLM
    for iter=1:maxiter
        %sumo las componentes del vector mediante el comando sum
        suma_y=ytotal;
        suma_alfa=alfatotal; 
        %Calculamos las matrices a1 y a2
        a1=n*suma_alfa^2*suma_y.^(n-1);
        a2=2*suma_alfa*suma_y.^n;
        %Matrices de diferenciacion
        dy=D1*suma_y;
        ddy=D2*suma_y;
        %Calculo el residuo
        r=-(ddy+(2./t).*dy+suma_alfa^2*suma_y.^n);
        %Contruimos la matriz A y B
        A=D2+diag(2./t)*D1+diag(a1);
        B=a2;
        %Elimino primera y última fila/columna
        Amatriz=A(2:N,2:N);
        Bmatriz=B(2:N);
        Rmatriz=r(2:N);
        %Defino las matrices aumentadas
        Asistema=[Amatriz,Bmatriz;D(N+1,2:N),0];
        Rsistema=[Rmatriz;0];
        %Resuelvo el sistema aumentado de todas las matrices juntas
        solucion=Asistema\Rsistema;
        %Saco la solucion una a una del sistema aumentado
        Y=zeros(N+1,1);
        Y(2:N)=solucion(1:end-1);
        alfa_i=solucion(end);
        %Actualizo los vectores
        ytotal=ytotal+Y;
        alfatotal=alfatotal+alfa_i;
        alfas(iter)=sum(alfatotal);
    end
    %Soluciones finales haciendo el sumatorio
    alfa=alfatotal;
    y=ytotal;
end
function D=chebdiff(N)
%CHEBDIFF Matriz de diferenciación de Chebyshev.
%D=CHEBDIFF(N) calcula la matriz de diferenciación de Chebyshev de tamaño (N+1)x(N+1).
%Puntos de Chebyshev-Gauss-Lobatto
z=cos(pi*(0:N)'/N);
%Coeficientes para los polinomios de Chebyshev
c=[2;ones(N-1,1);2];
c=c.*(-1).^(0:N)';
%Matriz de diferenciación, elementos D(j,K)
X=repmat(z,1,N+1);
dX=X-X';
D=(c*(1./c)')./dX;
D(eye(N+1)~=0)=0; %Se anulan los elementos diagonales
%Nuevos elementos diagonales
for k=1:N-1
    D(k+1,k+1)=-z(k+1)/(2*(1-z(k+1)^2));
end
%Condiciones para 00 y NN
D(1,1)=(2*N^2+1)/6;
D(N+1,N+1)=-D(1,1);
end
  