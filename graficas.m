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

% Este archivo es para dibujar gráficas.
N = 40; % Número de puntos para la malla.
maxiter = 10; % Iteraciones del método.
n_valores = [0 1 2 3 4]; % Valores de n.
% Figura 1.
% Gráfica comparativa de los métodos SLM y QLM.
figure % Para que nos muestre la gráfica. 
hold on 
grid on 
legend_entries = {}; % Esto me lo ha dicho ChatGPT, es para almacenar la leyenda.
for n = n_valores % Bucle para los valores de n.
    % Método QLM
    [alfa_q, y_q, ~, ~, x_q] = laneemdenNEWTONKANTOROVICH(n,N,maxiter);
    xi_q = alfa_q * x_q; % Para dibujarlo en la variable adimensional.
    % Método SLM 
    [alfa_s, y_s, ~] = laneemdenSLM(n,N,maxiter);
    z = cos(pi*(0:N)'/N); % Puntos de la malla.
    x_s = (z+1)/2; % Para convertir la malla a la misma que la del QLM.
    xi_s = alfa_s * x_s; % Para dibujarlo en la variable adimensional.
    % Colores para los distintos n.
    switch n % Para poner un color a cada n.
        case 0, col = 'k';
        case 1, col = 'b';
        case 2, col = [0 0.6 0];
        case 3, col = 'm';
        case 4, col = [0.85 0.33 0.1];
    end
    % Gráfica para el método QLM.
    plot(xi_q, y_q, '--', 'Color', col, 'LineWidth', 1.8)
    legend_entries{end+1} = sprintf('QLM n=%d', n);
    % Gráfica para el método SLM.
    plot(xi_s, y_s, ':', 'Color', col, 'LineWidth', 2)
    legend_entries{end+1} = sprintf('SLM n=%d', n);
end
% Nombre de los ejes y título.
xlabel('\xi')
ylabel('\theta(\xi)')
title('Ecuación de Lane-Emden: soluciones numéricas')
% Para poner los ejes entre los puntos que se quieran.
xlim([0 20])
ylim([-0.2 1.05])
% Leyenda de la gráfica.
legend(legend_entries,'Location','northeastoutside')

% Gráfica de las soluciones analíticas de la ecuación de Lane-Emden.
% Figura 2.
figure
xi = linspace(0,20,1000); % Eje x equiespaciado.
% Soluciones para los n conocidos.
theta_n0 = 1 - xi.^2/6;
theta_n1 = sin(xi)./xi;
theta_n1(1) = 1; % Para evitar dividir entre cero.
theta_n5 = (1 + xi.^2/3).^(-1/2);
% Gráfica.
hold on
plot(xi, theta_n0, 'LineWidth',2)
plot(xi, theta_n1, 'LineWidth',2)
plot(xi, theta_n5, 'LineWidth',2)
hold off
% Para poner los ejes entre los puntos deseados. 
xlim([0 20])
ylim([-2 1.05])
% Leyenda, título y nombre de los ejes.
legend('n=0','n=1','n=5','Location', 'southeast')
xlabel('\xi') 
ylabel('\theta(\xi)') 
title('Soluciones explícitas de Lane–Emden')
grid on

% Gráfica comparativa entre RK4 y ode45 para Lane-Emden.
% Figura 3.
epsilon = 1e-6;% Para evitar la singularidad en xi = 0.
puntofinal = 20;% Dominio en xi.
h_rk = 0.5;% Paso fijo RK4 que correspode con N=40.
% Muestra la figura.
figure
hold on
grid on
legend_entries = {};
for n = n_valores
    % Método RK4.
    [xi_rk, theta_rk] = laneemdenRK4(n, epsilon, puntofinal, h_rk);
    % Función ode45.
    [xi_ode, theta_ode] = laneemdenode45(n, epsilon, puntofinal);
    % Colores para cada n.
    switch n
        case 0, col = 'k';
        case 1, col = 'b';
        case 2, col = [0 0.6 0];
        case 3, col = 'm';
        case 4, col = [0.85 0.33 0.1];
    end
    % Gráfica RK4.
    plot(xi_rk, theta_rk, '--', 'Color', col, 'LineWidth', 1.6)
    legend_entries{end+1} = sprintf('RK4 n=%d', n);
    % Gráfica ode45.
    plot(xi_ode, theta_ode, '-', 'Color', col, 'LineWidth', 2)
    legend_entries{end+1} = sprintf('ode45 n=%d', n);
end
% Nombre de los ejes y título.
xlabel('\xi')
ylabel('\theta(\xi)')
title('Ecuación de Lane-Emden: RK4 vs ode45')
% Ajuste de los ejes.
xlim([0 20])
ylim([-0.2 1.05])
% Leyenda.
legend(legend_entries, 'Location', 'northeastoutside')

