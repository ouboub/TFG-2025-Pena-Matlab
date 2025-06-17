function [alfa,a,x,y]=laneemdenNEWTONKANTOROVICH(n,N,maxiter)
    % Resuelve la ecuación de Lane-Emden usando método espectral de Chebyshev
    % y el de Newton-Kantorovich
    % ENTRADAS:
    %   n: exponente en la ecuación Lane-Emden
    %   N: número de puntos de colocación
    %   maxiter: máximo número de iteraciones de Newton
    % SALIDAS:
    %   alfa: autovalor (solución de la ecuación)
    %   a: coeficientes de Chebyshev
    %   x: puntos de colocación
    %   y: solución evaluada en x
    
    %1.Puntos de Chebyshev-Lobatto en [0,1]
    j=(0:N-1)';
    x=cos(pi*j/(N-1));  %Puntos en [-1,1]
    x=(1+x)/2;        %Transformacion a [0,1]
    x(1)=1; 
    x(end)=0; %Asegurar exactamente 1 y 0 (evitar singularidades)    
    %2.Matrices de diferenciación de Chebyshev
    
    
    %3.Condiciones iniciales
    y0=cos((pi/2)*x); 
    alfa0=3;
    %4.Encontrar coeficientes iniciales de Chebyshev
   
  
    %5.Iteraciones de Newton-Kantorovich
%aqui iria el bucle
    %6.Evaluar la solucion final
end