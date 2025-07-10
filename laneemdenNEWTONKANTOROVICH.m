function [alfa,y,alfas]=laneemdenNEWTONKANTOROVICH(n,N,maxiter)
    %   Resuelve la ecuación de Lane-Emden usando método espectral de Chebyshev
    %   y el de Newton-Kantorovich
    %   ENTRADAS:
    %   n: exponente en la ecuación Lane-Emden % UB:10.07.2025:08:17: OK [0 5]
    %   N: número de puntos de colocación % UB:10.07.2025:08:17:OK
    %   maxiter: máximo número de iteraciones de Newton % UB:10.07.2025:08:17: OK
    %   SALIDAS:
    %   alfa: autovalor (solución de la ecuación) % UB:10.07.2025:08:16: OK
    %   y: solución evaluada en x% UB:10.07.2025:08:17: no entiendo
    %   alfas: distintas soluciones para cada n % UB:10.07.2025:08:17:no entiendo
    %   Falta o no está explicada (solo para saber donde esta su fallo)
    %   Y:            Valor a la malla 
    %   XCheb:    Puntos de la malla
    %   a:        Los coefficientes de la seria de Chebyshev

    


    % 1.Puntos de Chebyshev-Lobatto en [0,1]
    j=(0:N-1)';
    tj=pi*(j/(N-1));  % Puntos en [-1,1]
    x=(1+cos(tj))/2;  % Transformacion a [0,1]    
    % 2.Matrices de diferenciación de Chebyshev
    [D0,D1,D2]=chebdiff(N); % UB:10.07.2025:08:21: que es esto? estas matrices debe calcular usted
    % 3.Condiciones iniciales
    y0=cos((pi/2).*x); 
    alfa0=3;
    % 4.Encontrar coeficientes iniciales de Chebyshev
    a=D0\y0;
    % 5.Iteraciones de Newton-Kantorovich
    alfas=zeros(maxiter,1);
    alfas(1)=alfa0;
    alfa=alfa0;
    for iter=1:maxiter
        % Evaluaamos y y sus derivadas
        y=D0*a;
        yx=D1*a;
        yxx=D2*a;
        % Calculamos el residuo
        r=yxx+(2./x).*yx+alfa^2*y.^n; % puede haber fallo 
        % Iniciamos el sistema lineal para la corrección (Jacobiano)
        J=zeros(N+1,N+1);
        R=zeros(N+1,1);
        % Ecuación en puntos interiores (N-2 ecuaciones)
        for i=2:N-1
            J(i-1,1:N)=D2(i,:)+(2/x(i))*D1(i,:)+alfa^2*n*y(i)^(n-1)*D0(i,:);
            J(i-1,N+1)=2*alfa^n*y(i); % puede haber fallo
            R(i-1)=-r(i);
        end
        % Condiciones de frontera
        % y(0)=1
        J(N-1,1:N)=D0(1,:);
        R(N+1)=-(y(1)-1);
        % yx(0)=0 
        J(N,1:N)=D1(end,:);
        % y(1)=0
        J(N+1,1:N)=D0(end,:);
        % Resolver el sistema lineal
        delta=J\R;
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
    tj=pi*((0:N-1)'/(N-1));
    % Construimos D0 usando la relación trigonométrica T_n(cos(t))=cos(nt)
    D0=zeros(N,N);
    for i=1:N
        t=tj(i);
        for j=0:N-1
            D0(i,j+1)=cos(j*t);  % Equivalente a T_n(2x-1)
        end
    end
    % Construimos D1 y D2 usando fórmulas
    D1=zeros(N,N);
    D2=zeros(N,N);
    for i=1:N
        t=tj(i);
        for j=0:N-1
            % Primera derivada
            D1(i,j+1)=2*j*sin(j*t)/sin(t);
            D1(1,j+1)=2*j^2;
            D1(N,j+1)=2*(-1)^j*j^2;
            % Segunda derivada
            D2(i,j+1)=4*(-j^2*cos(j*t)/(sin(t)^2)+j*sin(j*t)*cos(t)/(sin(t)^3));
            D2(1,j+1)=(4/3)*j^2*(j^2-1);
            D2(N,j+1)=(4)*(-1)^j*(j^2)*((j^2-1)/3);
        end
    end
end
