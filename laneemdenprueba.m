%En este archivo vamos a crear una tabla para expresar los resultados
%obtenidos de los métodos QLM y SLM, para los distintos valores de n y la
%iteración en la que nos encontramos
N=40; %número de puntos para la malla
maxiter=10; %Número de iteraciones
nvalores=[0,1,2,3,4,5];  %Valores del índice politrópico
%Prepara la tabla de resultados
resultado=zeros(maxiter,length(nvalores));
    %Calcula resultados para cada n
    for i=1:length(nvalores)
        n=nvalores(i);
        %pongo el método que quiera probar
        [~,~,alfas]=laneemdenSLM(n,N,maxiter);
        resultado(:,i)=alfas;
    end
%Muestra la tabla formateada
disp('-------------------------------------------------------------');
disp('| Iters | n = 0       | n = 1       | n = 2       | n = 3       | n = 4    | n = 5       |');
disp('-------------------------------------------------------------');
    for iter=1:maxiter
        fprintf('| %-8d ',iter);
        for columnas=1:length(nvalores)
            fprintf('| %-20.15f ',resultado(iter,columnas));
        end
        fprintf('|\n');
    end
disp('-------------------------------------------------------------');

