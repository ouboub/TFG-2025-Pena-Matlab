% Parámetros de la simulación
% N=input('introduce N');    %Número de puntos de colocación
% maxiter=input('introduce maxiter');   %Número de iteraciones
N=40;
maxiter=10;
nvalores=[0,1,2,3,4,5];  %Valores del índice politrópico
%Prepara la tabla de resultados
resultado=zeros(maxiter,length(nvalores));
    %Calcula resultados para cada n
    for i=1:length(nvalores)
        n=nvalores(i);
        [~,~,alfas]=laneemdenSL(n,N,maxiter);
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

