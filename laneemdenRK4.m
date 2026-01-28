function [t,theta]=laneemdenRK4(n,epsilon,puntofinal,N)
    % Método RK4 para resolver la ecuacion de Lane-Emden
    % ENTRADAS:
    %   n: Índice politrópico
    %   epsilon: Valor inicial pequeño
    %   puntofinal: punto final de integración
    %   N: Paso fijo para RK4
    % SALIDAS:
    %   t: punto donde se calcula la solucion.
    %   theta: solucion de la ecuacion.

    %Inicialización
    t=epsilon:N:puntofinal; %vector del paso
    M=length(t);
    y=zeros(M,2); %matriz de la solucion
    theta0=1;
    dtheta0=0;
    y(1,:)=[theta0,dtheta0]; %condiciones iniciales
    funcion=@(t,y)[y(2);-(2/t)*y(2)-y(1)^n]; %ecuacion de lane-emden
    %Iteración RK4
    for i=1:M-1
        ti=t(i);
        yi=y(i,:)';
        %Etapas RK4
        k1=funcion(ti,yi);
        k2=funcion(ti+N/2,yi+(N/2)*k1);
        k3=funcion(ti+N/2,yi+(N/2)*k2);
        k4=funcion(ti+N,yi+N*k3);
        %Actualización
        y(i+1,:)=yi+(N/6)*(k1+2*k2+2*k3+k4);
    end
%Saco la soluciones
theta=y(:,1);
%Si se alcanza theta=0 (es decir la solucion) antes del puntofinal
%la solución la mostrará haciendo lo siguiente, que es cuando cruza el primer cero 
%aquí me he ayudado de chatgpt porque no sabia como podia programar que la
%me diera la solucion cuando cruza
zero_crossing=find(theta <= 0, 1); 
if ~isempty(zero_crossing)           
    t=t(zero_crossing(end));
    theta=theta(zero_crossing(end));
end
%para poder dibujar la grafica necesitariamos comentar estas últimas cinco
%líneas ya que si no solo dibujaria un punto, el de la solución