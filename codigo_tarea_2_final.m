%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Evolución de H_rms con roller (Lippmann et al.,1996) y sin roller, pero con rotura
%  (Alsina & Baldock, 2007)
%
%
%  //////////////////////////// CON ROLLER ///////////////////////////////
%  Ecuaciones del modelo que estoy usando:
%    dF/dx = - Qb * eps_r                         
%    F     =  Ew * Cg  +  Qb * Er                 
%
%  Donde:
%    Ew    = (rho*g/8) * Hrms^2                   (energía de ola)
%    Er    = (1/8) * rho * c * f * Hb^3/(h*tanσ)  (energía del roller)
%    eps_r = (1/4) * rho * g * f * Hb^3/h * cosσ  (disipación del roller)
%    Qb    = exp( - (Hb/Hrms)^2 )                 (prob. rompimiento Rayleigh)
%
%  Hb, sacado de Alsina 2007:
%    Hb = (0.88/k) * tanh( (gamma_b * k * h) / 0.88 )
%
%  Otros:
%    omega^2 = g*k*tanh(k*h),  c = omega/k,
%    Cg = 0.5*c*(1 + 2kh/sinh(2kh))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------Comentarios sobre ROLLER---------------------: 
%Me di cuenta que si utilizo Hb^2 en Er y epsilonr, la línea de la altura de
%ola se ve más parabólica y no tan "plana" como la que entrega este código, ahí calzaría
%mejor con los puntos R37 y R39, pero según Lippmann es Hb^3. Si no me equivoco Baldock usa Hb^2.
%En comparación con la tarea 1, la Hrms baja al acercarse a la orilla gracias a la disipación roller. 
%
%La idea es ir variando gamma_b y sigma_deg. A menor gamma ~ 0, la curva de
%Hrms ya llegando a la orilla tiende a levantarse y se aleja de lo
%esperado. Con gamma ~ 1, Hrms baja llegando a la orilla.
%
%Al menos para R39 (parámetros T=8 y Hrms0=0.37) se ven buenos resultados
%para los primeros 3 puntos, sin embargo no logré hacer que la curva bajara
%para calzar con los 3 restantes. Hice otro código (no adjunto) con Hb^2 (Baldock) y
%ahí calza bien por su forma parabólica, pero quise seguir con lo que propone Lippmann. 

clear; 
clc;    

%% ---------------------- Condiciones de ola (RXX) ---------------------------
archivo   = 'batimetria.txt';   
%T y Hrms0 variar de acuerdo a experimento elegido. En este caso se varía
%entre T = 8.0 y 5.0, y Hrms0 = 0.37 y 0.51 
T         = 8.0;              %R39  
Hrms0     = 0.37;             %R39  
g         = 9.81;              
rho       = 1025;               

%% ---------------------- Parámetros del modelo --------------------------
gamma_b     = 0.99;       %probar valores para ajustar     
sigma_deg   = 2.0;        %probar valores para ajustar   

%%%% ---------------------- Otros para Newton --------------------------
h_min_phys  = 1e-3;            
tol         = 1e-12;           
maxit       = 60;              
x_stop = 86.62;                

%% ---------------------- Lectura de la batimetría ------------------------------
% el archivo batimetria.txt tiene dos columnas: x y z.
% por convención z<0 bajo el nivel medio del mar (profundidad positiva h=-z).
datos = readmatrix(archivo);                 
if size(datos,2) < 2
    error('El archivo debe tener dos columnas: x y z.');
end
x = datos(:,1);                              
z = datos(:,2);                              

N  = numel(x);                              
dx = diff(x); dx(end+1) = dx(end);          
h  = max(0, -z);                             

i_stop = find(x >= x_stop, 1, 'first');
if isempty(i_stop)
    i_stop = N;  
end

%% ---------------------- Cinemática de la ola ------------------------------
omega = 2*pi/T;     
f     = 1/T;         

k  = zeros(N,1);    
C  = zeros(N,1);     
Cg = zeros(N,1);     

for i = 1:N
    hi = max(h(i), h_min_phys);  
    
    % semilla inicial para k. Uso la de aguas profundas o medias.
    k_i = max(omega^2/g, omega/sqrt(g*hi + eps));
    
    % método de Newton para resolver omega^2 = g*k*tanh(k*h)
    for it = 1:maxit
        th    = tanh(k_i*hi);          
        sech2 = 1/cosh(k_i*hi);        
        sech2 = sech2*sech2;           
       
        Fk  = omega^2 - g*k_i*th;                  
        dFk = -g*( th + k_i*hi*sech2 );            
        k_new = k_i - Fk/dFk;                     
        
        if abs(k_new - k_i) < tol
            k_i = k_new;                          
            break
        end
        k_i = k_new;                               
    end
    
    k(i)  = k_i;                  
    C(i)  = omega/k_i;            
    Cg(i) = 0.5*C(i)*(1 + (2*k_i*hi)/sinh(2*k_i*hi));  
end

%% ---------------------- Roller -----------------
% la altura máxima de ola antes de romper (aparece en Alsina)
Hb = (0.88./k) .* tanh( (gamma_b .* k .* h) ./ 0.88 );

% se pasa el sigma_deg a radianes
tan_sig = tan(deg2rad(sigma_deg));   
cos_sig = cos(deg2rad(sigma_deg));   

% energía del roller Er y disipación eps_r
% se usa max() para evitar dividir por cero cuando h tiene a 0
Er    = (1/8) * rho .* C .* f .* (Hb.^3) ./ max(h.*tan_sig, 1e-12);   
eps_r = (1/4) * rho * g * f .* (Hb.^3) ./ max(h,1e-12) .* cos_sig;    

%% ---------------------- PROBABILIDAD DE ROMPIMIENTO Qb ------------------
% Qb = exp( - (Hb/Hrms)^2 ), Alsina también
Qb      = @(Hr,Hb_) exp( - (Hb_ ./ max(Hr,1e-12)).^2 );
dQb_dHr = @(Hr,Hb_) Qb(Hr,Hb_) .* ( 2*(Hb_.^2) ./ max(Hr,1e-12).^3 );

%% ---------------------- INTEGRACIÓN DEL BALANCE DE ENERGÍA -------------
Hrms = NaN(N,1);     
F    = NaN(N,1);    

% condición inicial
Hrms(1) = Hrms0;                       
Ew1     = (rho*g/8) * Hrms(1)^2;       
Qb1     = Qb(Hrms(1), Hb(1));          
F(1)    = Ew1*Cg(1) + Qb1*Er(1);       

for i = 1 : min(N-1, i_stop-1)
    
    % si prácticamente no hay agua, fijo cosas para no explotar
    if h(i) < h_min_phys
        F(i+1)    = max(F(i), 0);      
        Hrms(i+1) = 0;                
        continue
    end
    
    % disipación Dr = Qb * eps_r
    Qbi = Qb(Hrms(i), Hb(i));          
    Dr  = Qbi * eps_r(i);          
    
    % Avanzo el flujo con un paso explícito (dF/dx = ΔF/Δx)
    F(i+1) = max(F(i) - Dr*dx(i), 0); 
    
    % Hrms(i+1)
    % F_{i+1} = (rho g/8) H^2 * Cg(i+1) + Qb(H, Hb(i+1)) * Er(i+1)
    cgn1 = Cg(i+1);
    Ern1 = Er(i+1);
    Hb1  = Hb(i+1);
    
    % parto desde el valor anterior
    H = max(Hrms(i), 1e-6);
    
    for it = 1:maxit
        QbH = Qb(H, Hb1);                                         
        G   = (rho*g/8)*H^2*cgn1 + QbH*Ern1 - F(i+1);            
        if abs(G) < 1e-10
            break                                                
        end
        dQ  = dQb_dHr(H, Hb1);                                    
        dG  = (rho*g/4)*H*cgn1 + dQ*Ern1;                        
        dG  = sign(dG)*max(abs(dG),1e-12);                        
        Hn  = max(H - G/dG, 0);                                  
        if abs(Hn - H) < 1e-10
            H = Hn; break
        end
        H = Hn;                                                   
    end
    Hrms(i+1) = H;                                               
end

%% ---------------------- GRÁFICO DE RESULTADOS --------------------------
figure('Color','w'); hold on; grid on;

% Fondo (z) en azul
hBathy = plot(x, z, 'b-', 'LineWidth', 2);

% Curva de Hrms en rojo (NaN después del corte no se dibuja)
hHrms  = plot(x, Hrms, 'r-', 'LineWidth', 2);

% Línea del nivel 0 (superficie media)
yline(0,'k-');

xlabel('x [m]');
ylabel('Cota z y H_{rms} [m]');
title(sprintf(['Evolución H_{rms} con roller (Lippmann)  ', ...
               'T=%.2fs, H_{rms,0}=%.2fm, \\gamma=%.2f, \\sigma=%g^\\circ'], ...
               T, Hrms0, gamma_b, sigma_deg));

% --- Puntos de comparación (R37 en rosado, R39 en verde) ---
x_R37 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R37 = [0.496623780748516 0.510497239987004 0.497355602540387 ...
         0.269045854415693 0.267893303559896 0.212864134792817];
hR37 = scatter(x_R37, H_R37, 70, 'o', ...
    'MarkerFaceColor', [1 0.2 0.6], 'MarkerEdgeColor', 'k');

x_R39 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R39 = [0.434633492839042 0.444342173431960 0.439117550441414 ...
         0.273265428250874 0.276439604842855 0.211498629675050];
hR39 = scatter(x_R39, H_R39, 70, 's', ...
    'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', 'k');

legend([hBathy, hHrms, hR37, hR39], ...
       {'Cota z','H_{rms}','R37','R39'}, 'Location','best');

% Limito el eje X hasta el punto de corte para que el gráfico termine ahí
xlim([x(1) x(i_stop)]);

% Ajusto el eje Y para que se vea todo razonable
ymin = min([z; -0.05]);                      % bajo un poco de 0 para ver el eje
ymax = max([z(~isnan(z)); Hrms(~isnan(Hrms))]) * 1.15;
ylim([ymin, ymax]);

%% =======================================================================
%
% ////////////////////////// CON ROTURA Y SIN ROLLER //////////////////////
%
% Alsina & Baldock (2007)
%
% d/dx(E*Cg) = -D, con D según ecuación analítica (Eq. 11 del paper).
% batimetria.txt: [x(m), z(m)] con z NEGATIVO bajo el mar y POSITIVO en seco.
clear; clc;

% ------------------ Entradas ------------------
archivo = 'batimetria.txt';
T      = 8.0;                 
Hrms0  = 0.37;                
g      = 9.81;                
rho    = 1025;                
theta  = 0;                   
tol    = 1e-12;              
maxit  = 50;                  

% ---------- Parámetros (variar!) ---------------
B            = 1.0;           % parámetro de ajuste en A
gamma_fixed  = 0.8;          % índice de rompiente fijo
h_min_phys   = 1e-3;          % profundidad mínima física para el modelo [m]

% ------------------ Lectura -------------------
datos = readmatrix(archivo);
if size(datos,2) < 2
    error('El archivo debe tener dos columnas: x y z (cota).');
end
x = datos(:,1);
z = datos(:,2);               

N  = numel(x);
dx = diff(x); dx(end+1) = dx(end);

% profundidad para la ola (>=0). Si z>0 (seca), h=0.
h = max(0, -z);

% ------------------ PRE-CÁLCULOS DE OLA ------------------
omega = 2*pi/T;
fp    = 1/T;

k  = zeros(N,1);
C  = zeros(N,1);
Cg = zeros(N,1);

% resolver k(x) con Newton–Raphson (dispersión lineal)
for i = 1:N
    hi = max(h(i), h_min_phys);   % evita singularidades
    % Semilla
    if omega^2*hi/g > pi
        k_i = omega^2/g;
    else
        k_i = omega/sqrt(g*hi + eps);
    end
    % Newton
    for it = 1:maxit
        th   = tanh(k_i*hi);
        sech = 1/cosh(k_i*hi);
        F  = omega^2 - g*k_i*th;
        dF = -g*( th + k_i*hi*(sech^2) );
        k_new = k_i - F/dF;
        if abs(k_new - k_i) < tol, k_i = k_new; break; end
        k_i = k_new;
    end
    k(i)  = k_i;
    C(i)  = omega / k_i;
    Cg(i) = 0.5*C(i)*(1 + (2*k_i*hi)/sinh(2*k_i*hi)) * cos(theta);
end

% ---- Hb (umbral de rompiente) ----
gamma = gamma_fixed * ones(N,1);
Hb = (0.88./k) .* tanh( (gamma .* k .* max(h, h_min_phys)) ./ 0.88 );

% ------------------ INTEGRACIÓN DE FLUJO (dF/dx = -D) -------------------
A = 0.25 * rho * g * fp * B;           % (Ecuación 10)
Hrms = zeros(N,1); Hrms(1) = Hrms0;

F = zeros(N,1);
F(1) = (rho*g/8) * Hrms(1)^2 * Cg(1);

for i = 1:N-1
    % Si estamos ya prácticamente en seco, detenemos propagación
    if h(i) < h_min_phys
        F(i+1)    = max(F(i), eps);
        Hrms(i+1) = 0;
        continue
    end

    % ---- NUEVA D (Ecuación 11) evaluada en el nodo i ----
    q = Hb(i) / max(Hrms(i), 1e-12);
    bracket = (q^3 + 1.5*q) * exp(-(q^2)) + (3/4)*sqrt(pi)*(1 - erf(q));
    Di = A * (Hrms(i)^3) / max(h(i),1e-12) * bracket;   % [W/m^2]
    
    % Avance de flujo y actualización de Hrms
    F(i+1) = max(F(i) - Di*dx(i), eps);
    Hrms(i+1) = sqrt( max( 8*F(i+1) / (rho*g*max(Cg(i+1),eps)) , 0 ) );
end

% ------------------ GRÁFICO ÚNICO ------------------
figure('Color','w'); hold on; grid on;

% Línea de cota (z) -> positiva en seco, negativa bajo el mar
hBathy = plot(x, z, 'b-', 'LineWidth', 2);      % ahora graficamos z(x) directamente
% Hrms (siempre positivo)
hHrms  = plot(x, Hrms, 'r-', 'LineWidth', 2);
yline(0,'k-');                                   % nivel medio del mar

% ---- Puntos R37 (ROSADOS) ----
x_R37 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R37 = [0.496623780748516 0.510497239987004 0.497355602540387 ...
         0.269045854415693 0.267893303559896 0.212864134792817];
hR37 = scatter(x_R37, H_R37, 70, 'o', ...
    'MarkerFaceColor', [1 0.2 0.6], 'MarkerEdgeColor', 'k');

% ---- Puntos R39 (VERDES) ----
x_R39 = [23.45 45.45 52.70 60.20 70.95 81.95];
H_R39 = [0.434633492839042 0.444342173431960 0.439117550441414 ...
         0.273265428250874 0.276439604842855 0.211498629675050];
hR39 = scatter(x_R39, H_R39, 70, 's', ...
    'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', 'k');

xlabel('x [m]');
ylabel('Cota z y H_{rms} [m]');
title(sprintf('Evolución H_{rms} con rotura (Alsina & Baldock)  T=%.2fs, H_{rms,0}=%.2fm', T, Hrms0));
legend([hBathy, hHrms, hR37, hR39], {'Cota z','H_{rms}','R37','R39'}, 'Location','best');

% Límites: cubre el mínimo de z y el máximo entre z y Hrms
ymin = min([z; -0.05]);
ymax = max([z; Hrms]) * 1.15;
ylim([ymin, ymax]);





