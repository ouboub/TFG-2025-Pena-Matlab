N = 40; %puntos para la malla
maxiter = 10; %iteraciones del metodo
n_valores = [0 1 2 3 4]; %valores de n
figure
hold on
grid on
legend_entries = {}; %esto me lo ha dicho chatgpt, es para almacenar la leyenda
for n = n_valores
    %método QLM
    [alfa_q, y_q, ~, ~, x_q] = laneemdenNEWTONKANTOROVICH(n,N,maxiter);
    xi_q = alfa_q * x_q; %para conseguir que esté en la variable adimensional
    %método SLM 
    [alfa_s, y_s, ~] = laneemdenSLM(n,N,maxiter);
    z = cos(pi*(0:N)'/N);
    x_s = (z+1)/2; %para convertir la malla a la misma que la del QLM
    xi_s = alfa_s * x_s; %para conseguir que esté en la variable adimensional
    %colores para los distintos n
    switch n
        case 0, col = 'k';
        case 1, col = 'b';
        case 2, col = [0 0.6 0];
        case 3, col = 'm';
        case 4, col = [0.85 0.33 0.1];
    end
    %Gráficas de los dos metodos
    plot(xi_q, y_q, '--', 'Color', col, 'LineWidth', 1.8)
    legend_entries{end+1} = sprintf('QLM n=%d', n);

    plot(xi_s, y_s, ':', 'Color', col, 'LineWidth', 2)
    legend_entries{end+1} = sprintf('SLM n=%d', n);

end
xlabel('\xi')
ylabel('\theta(\xi)')
title('Ecuación de Lane-Emden: soluciones numéricas')
%Para poner los ejes entre los puntos que se quieran
xlim([0 20])
ylim([-0.2 1.05])
%leyenda de la gráfica
legend(legend_entries,'Location','northeastoutside')
