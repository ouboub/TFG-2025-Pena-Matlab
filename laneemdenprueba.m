% Parámetros de la simulación
N=input('introduce N');    %Número de puntos de colocación
maxiter=input('introduce maxiter');   %Número de iteraciones
nvalores= [1, 2, 3, 4];  %Valores del índice politrópico
%Prepara la tabla de resultados
resultado=zeros(maxiter,length(nvalores));
    %Calcula resultados para cada n
    for i=1:length(nvalores)
        n=nvalores(i);
        [~, ~,alfas]=laneemden(n,N,maxiter);
        resultados(:,i)=alfas;
    end
%Muestra la tabla formateada
disp('-------------------------------------------------------------');
disp('| Iters | n = 1       | n = 2       | n = 3       | n = 4       |');
disp('-------------------------------------------------------------');
    for iter=1:maxiter
        fprintf('| %-5d ',iter);
        for columnas=1:length(nvalores)
            fprintf('| %-12.6f ',resultados(iter,columnas));
        end
        fprintf('|\n');
    end
disp('-------------------------------------------------------------');

