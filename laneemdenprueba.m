% UB:02.06.2025:17:53: por favor usar input!
% UB:02.06.2025:17:53: ese programa ha sido generado por ChatGTP???
% Parámetros de la simulación

%% eso es una celda
% N=input('introduce N');    %Número de puntos de colocación
% maxiter=input('introduce maxiter');   %Número de iteraciones
N=70;
maxiter=10;
nvalores=[1,2,3,4];  %Valores del índice politrópico
%Prepara la tabla de resultados
resultado=zeros(maxiter,length(nvalores));
    %Calcula resultados para cada n
    for i=1:length(nvalores)
        n=nvalores(i);
        [~,~,alfas]=laneemdenNEWTONKANTOROVICH(n,N,maxiter);
        resultado(:,i)=alfas;
    end
%Muestra la tabla formateada
disp('-------------------------------------------------------------');
disp('| Iters | n = 1       | n = 2       | n = 3       | n = 4       |');
disp('-------------------------------------------------------------');
    for iter=1:maxiter
        fprintf('| %-5d ',iter);
        for columnas=1:length(nvalores)
            fprintf('| %-12.6f ',resultado(iter,columnas));
        end
        fprintf('|\n');
    end
disp('-------------------------------------------------------------');

