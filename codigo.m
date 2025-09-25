% === Shoaling sin rotura sobre perfil con barra ===
% Lee batimetria.txt: [x(m) , h(m)] con h NEGATIVO bajo el mar
clear; clc;

% ------------------ ENTRADAS ------------------
archivo = 'batimetria.txt';
T = 5.0;                 % Periodo [s]
Hrms0 = 0.51;            % Hrms en el primer punto [m]
g = 9.81;                % [m/s^2]
rho = 1025;              % kg/m^3 (no afecta: se cancela)
theta = 0;               % incidencia normal
tol = 1e-12;             % tolerancia Newton
maxit = 50;              % iteraciones max

% ------------------ LECTURA -------------------
datos = readmatrix(archivo);
if size(datos,2) < 2
    error('El archivo debe tener dos columnas: x y h.');
end
x = datos(:,1);
h = abs(datos(:,2)); % profundidad positiva, aunque archivo tenga h negativa

% Ordenar x en caso necesario
if any(diff(x) <= 0)
    [x,ord] = sort(x);
    h = h(ord);
end

% ------------------ PRE-CÁLCULOS --------------
omega = 2*pi/T;
N = numel(x);

k  = zeros(N,1);
C  = zeros(N,1);
Cg = zeros(N,1);

% Resolver k(x) con Newton–Raphson
for i = 1:N
    hi = h(i);

    % Semilla según régimen
    if omega^2*hi/g > pi
        k_i = omega^2/g;                 % aprox aguas profundas
    else
        k_i = omega/sqrt(g*hi + eps);    % aprox aguas someras
    end

    for it = 1:maxit
        th   = tanh(k_i*hi);
        sech = 1/cosh(k_i*hi);
        F  = omega^2 - g*k_i*th;
        dF = -g*( th + k_i*hi*(sech^2) );
        k_new = k_i - F/dF;
        if abs(k_new - k_i) < tol
            k_i = k_new; break;
        end
        k_i = k_new;
    end
    k(i) = k_i;
    C(i) = omega / k_i;
    Cg(i)= 0.5*C(i)*(1 + (2*k_i*hi)/sinh(2*k_i*hi)) * cos(theta);
end

% ------------------ PROPAGACIÓN ---------------
Hrms = zeros(N,1);
Hrms(1) = Hrms0;

for i = 1:N-1
    Hrms(i+1) = Hrms(i) * sqrt( Cg(i) / max(Cg(i+1),eps) );
end

% ------------------ GRÁFICO ------------------
figure('Color','w'); hold on; grid on;

% Curvas base
hBathy = plot(x, -h, 'b-', 'LineWidth', 2);      % batimetría (negativa)
hHrms  = plot(x, Hrms, 'r-', 'LineWidth', 2);    % Hrms (sobre 0)
hSea   = yline(0,'k-'); %#ok<NASGU>              % nivel del mar (no va a la leyenda)

% ---- Puntos R37 (ROSADOS) ----
x_R37 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R37 = [0.496623780748516 0.510497239987004 0.497355602540387 ...
         0.269045854415693 0.267893303559896 0.212864134792817];
hR37 = scatter(x_R37, H_R37, 70, 'o', ...
    'MarkerFaceColor', [1 0.2 0.6], ...   % rosado
    'MarkerEdgeColor', 'k');

% ---- Puntos R39 (VERDES) ----
x_R39 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R39 = [0.434633492839042 0.444342173431960 0.439117550441414 ...
         0.273265428250874 0.276439604842855 0.211498629675050];
hR39 = scatter(x_R39, H_R39, 70, 's', ...
    'MarkerFaceColor', [0 0.7 0], ...     % verde
    'MarkerEdgeColor', 'k');

xlabel('x [m]');
ylabel('Cota / H_{rms} [m]');
title(sprintf('Batimetría y evolución de H_{rms} (T=%.2f s, H_{rms,0}=%.2f m)', T, Hrms0));

% ¡Ojo! Pasamos los handles explícitamente para que la leyenda no tome otras líneas
legend([hBathy, hHrms, hR37, hR39], ...
       {'Batimetría','H_{rms}','R37','R39'}, 'Location','best');

ylim([-max(h)*1.2 , max(Hrms)*1.2]);
