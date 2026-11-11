function cfg = config_experiment()
%CONFIG_EXPERIMENT  Parámetros centralizados del estudio comparativo de
%   algoritmos de reconstrucción EIT diferencial (MATLAB + EIDORS).
%
%   cfg = config_experiment() devuelve una estructura con TODOS los
%   parámetros modificables del estudio. Ningún otro fichero contiene
%   valores numéricos de configuración: para cambiar el experimento basta
%   con editar este archivo (o modificar la estructura antes de llamar a
%   main_EIT_comparison(cfg)).
%
%   Convenciones:
%     - Coordenadas de los fantomas normalizadas al radio del dominio (R = 1).
%     - Conductividades en S/m.
%     - Datos diferenciales: dV = V1 - V0 (V1 con inclusiones, V0 homogéneo).
%     - Todos los algoritmos reconstruyen dsigma = sigma - sigma_fondo.

cfg = struct();

%% ------------------------------------------------------------ General
cfg.name           = 'EIT_comparison';
cfg.run_mode       = 'full';        % 'pilot' = solo prueba piloto; 'full' = estudio completo
cfg.eidors_startup = '';            % ruta a eidors/startup.m ('' si EIDORS ya está en el path)
cfg.output_root    = fullfile(pwd, 'results');
cfg.verbose        = true;
cfg.eidors_log_level = 1;           % nivel de mensajes de EIDORS (1 = mínimo)

%% ------------------------------------------------------------ Modelos FEM
cfg.model.n_elec    = 16;
cfg.model.fwd_mesh  = 'e2c2';       % malla fina: genera las medidas (mk_common_model)
cfg.model.inv_mesh  = 'c2c2';       % malla gruesa: reconstrucción (evita el inverse crime)
cfg.model.sigma_bkg = 1;            % S/m
cfg.model.stim.inj     = [0 1];     % excitación adyacente
cfg.model.stim.meas    = [0 1];     % medida adyacente
cfg.model.stim.options = {'no_meas_current'};   % sin medidas sobre electrodos de corriente
cfg.model.stim.amplitude = 1;       % A
cfg.model.radius_tol       = 0.02;  % tolerancia relativa en el radio de ambos modelos
cfg.model.elec_angle_tol_deg = 2;   % tolerancia angular entre electrodos de ambos modelos

%% ------------------------------------------------------------ Escenarios de prueba
% Cada escenario es un conjunto de inclusiones circulares:
%   center = [x y] normalizado, radius normalizado, sigma en S/m.
sc = struct('name', {}, 'label', {}, 'inclusions', {});
sc(1).name = 'S1_central';
sc(1).label = 'Inclusión central conductora';
sc(1).inclusions = struct('center', {[0 0]}, 'radius', {0.2}, 'sigma', {2});

sc(2).name = 'S2_periferica';
sc(2).label = 'Inclusión periférica conductora';
sc(2).inclusions = struct('center', {[0.6 0]}, 'radius', {0.2}, 'sigma', {2});

sc(3).name = 'S3_dos_inclusiones';
sc(3).label = 'Dos inclusiones conductoras próximas';
sc(3).inclusions = struct('center', {[-0.3 0], [0.3 0]}, ...
                          'radius', {0.2, 0.2}, 'sigma', {2, 2});

% Medio heterogéneo: una inclusión conductora y otra aislante a la vez
% (contraste de signo opuesto, como tejidos de distinta conductividad).
sc(4).name = 'S4_heterogeneo';
sc(4).label = 'Inclusiones conductora y aislante';
sc(4).inclusions = struct('center', {[-0.4 0.3], [0.4 -0.3]}, ...
                          'radius', {0.2, 0.2}, 'sigma', {2, 0.5});
cfg.scenarios = sc;

%% ------------------------------------------------------------ Ruido
% SNR en dB (Inf = sin ruido). El ruido se añade por separado a V0 y V1;
% su nivel se define según cfg.noise.reference (ver más abajo).
% Nivel de ruido en % (0 = sin ruido): con reference = 'difference', el RMS
% del ruido resultante en dV es ese porcentaje del RMS(dV limpio).
cfg.noise.pct            = [0 1 5 15];
cfg.noise.snr_db         = -20 * log10(cfg.noise.pct / 100);   % 0 % -> Inf (40, 26.0, 16.5 dB)
cfg.noise.labels         = {'N0', 'N1', 'N2', 'N3'};
cfg.noise.n_realizations = 10;      % realizaciones por nivel con ruido
cfg.noise.seed_base      = 1000;    % semilla base del conjunto de PRUEBA
% Referencia del RMS que define sigma_ruido = RMS_ref / 10^(SNR/20):
%   'difference'  -> RMS(dV limpio)/sqrt(2): el SNR indicado es el de dV
%                    (opción elegida: con 'acquisition' el SNR efectivo en dV
%                    a 30 y 20 dB quedaba entre -19 y +4 dB)
%   'acquisition' -> RMS de cada adquisición limpia (definición inicial)
cfg.noise.reference      = 'difference';

%% ------------------------------------------------------------ Conjunto de validación
% Fantomas aleatorios, independientes de los de prueba, para el ajuste
% de hiperparámetros. No se usan nunca en la evaluación.
cfg.val.n_phantoms        = 24;
cfg.val.seed              = 777;     % semilla de la geometría de validación
cfg.val.radius_range      = [0.12 0.25];
cfg.val.max_extent        = 0.85;    % |centro| + radio <= max_extent
cfg.val.max_inclusions    = 2;
cfg.val.sigma_choices     = [2 0.5]; % mismas conductividades que la prueba
cfg.val.min_gap           = 0.05;    % separación mínima entre inclusiones
cfg.val.n_realizations    = 2;       % realizaciones por nivel con ruido
cfg.val.seed_base         = 500000;  % semilla base del ruido de validación (disjunta de la de prueba)

%% ------------------------------------------------------------ Evaluación
cfg.eval.npix        = 32;          % cuadrícula común 32 x 32
cfg.eval.mask_radius = 0.95;        % región circular válida (fracción del radio)

%% ------------------------------------------------------------ Error por iteración
% Algoritmos iterativos cuyo RIE se registra en cada iterado (evaluación final)
cfg.iter_tracking.algorithms = {'TV', 'GN'};

%% ------------------------------------------------------------ Tiempos
cfg.timing.n_warmup = 3;            % reconstrucciones de calentamiento no registradas

%% ------------------------------------------------------------ Algoritmos
cfg.algorithms    = {'BP', 'NOSER', 'TSVD', 'TV', 'GN', 'GREIT'};
cfg.alg_labels    = {'BP (Sheffield)', 'NOSER', 'TSVD', 'TV', 'GN + Tikhonov', 'GREIT'};
% Calibración global de amplitud (un escalar por algoritmo, fijado en
% validación) para los métodos cuya salida no está en unidades físicas.
cfg.calibration.algorithms = {'BP', 'GREIT'};

% --- BP (Sheffield Mk I mediante mk_common_gridmdl('backproj'))
cfg.hp.BP.normalize   = [0 1];      % datos diferenciales normalizados (Sheffield) o no
cfg.bp.search_orientation = true;   % alinear numeración de electrodos (16 rotaciones x 2 reflexiones)

% --- NOSER (inv_solve_diff_GN_one_step + prior_noser)
cfg.hp.NOSER.lambda   = logspace(-4, 0, 17);
cfg.hp.NOSER.exponent = [0.25 0.5 0.75 1.0];

% --- TSVD (inv_solve_TSVD)
% Parámetro nativo de inv_solve_TSVD en EIDORS 3.12: PORCENTAJE de s_1.
% Se conservan los s_i con s_i/s_1 >= umbral/100 (verificado en la piloto).
% Con 16 electrodos el jacobiano tiene rango 104 y s_104/s_1 ~ 1.5e-6.
cfg.hp.TSVD.threshold = logspace(-3, 1.5, 19);    % 0.001 % ... 31.6 %
cfg.tsvd.verify_interpretation = true;            % identifica el significado del parámetro

% --- TV (inv_solve_TV_pdipm + prior_TV)
cfg.hp.TV.lambda    = logspace(-7, -2, 11);
cfg.hp.TV.max_iter          = [3 4 5 10];                    % antes [5 10 15]: el óptimo era 10, el mínimo probado
cfg.tv.term_tolerance = 1e-3;

% --- GN iterativo + Tikhonov
cfg.gn.implementation = 'explicit'; % 'explicit' (fwd_solve + calc_jacobian de EIDORS) o 'inv_solve_gn'
% Escala de referencia: s_1(J)^2 ~ 6e-3; lambda^2 debe recorrer el espectro de J'J.
cfg.hp.GN.lambda    = logspace(-5, -1, 17);
cfg.hp.GN.max_iter = [1 2 3 5 10];                   % antes [3 10]: el óptimo era 3, el mínimo probado
cfg.gn.tol_step     = 1e-3;         % parada: ||paso|| / ||dsigma|| < tol_step
cfg.gn.tol_cost     = 1e-4;         % parada: descenso relativo del coste < tol_cost
cfg.gn.max_backtrack = 6;           % reducciones del paso (búsqueda lineal)
cfg.gn.min_sigma    = 1e-3;         % conductividad mínima admisible (S/m)

% --- GREIT (mk_GREIT_model)
cfg.hp.GREIT.noise_figure = [0.5 1.0 2.0 3.0 4.0];   % fuera 0,3 (inalcanzable);
cfg.hp.GREIT.target_size    = [0.05 0.08 0.12 0.16];         % antes [0.05 0.08 0.12]
cfg.greit.radius    = 0.2;          % radio de la imagen deseada (fracción del dominio)
cfg.greit.Nsim      = 1000;         % objetivos de entrenamiento
cfg.greit.distr     = [];           % [] = distribución por defecto de EIDORS
cfg.greit.seed      = 4242;         % semilla del entrenamiento
cfg.greit.normalize = 0;            % mismos datos diferenciales que el resto
% mk_GREIT_model exige un modelo 3D: se entrena sobre la malla gruesa
% extruida a una losa de una capa con electrodos que ocupan todo el espesor.
cfg.greit.extrude_height = 0.04;    % espesor de la losa (fracción del radio)
cfg.greit.equiv_tol      = 0.05;    % error relativo máximo frente al modelo 2D
% Reintentos si mk_GREIT_model no alcanza el NF pedido con su tolerancia interna
cfg.greit.nf_retry_factors = [1.001 0.999 1.003 0.997 1.01 0.99];

%% ------------------------------------------------------------ Punto de control del ajuste
% Tras ajustar cada algoritmo se guarda el progreso; si la ejecución se
% interrumpe, la siguiente reanuda desde ahí (solo si rejillas, opciones y
% datos de validación son idénticos).
cfg.tuning.checkpoint_file = fullfile(cfg.output_root, 'tuning_checkpoint.mat');
cfg.tuning.resume          = true;

%% ------------------------------------------------------------ Prueba piloto
cfg.pilot.scenario      = 1;        % inclusión central, sin ruido
cfg.pilot.stop_on_fail  = true;
cfg.pilot.hp.BP    = struct('normalize', 1, 'transform', [0 1], 'scale', 1);
cfg.pilot.hp.NOSER = struct('lambda', 0.03, 'exponent', 0.5);
cfg.pilot.hp.TSVD  = struct('threshold', 0.1);   % valor por defecto de EIDORS (0,1 %)
cfg.pilot.hp.TV    = struct('lambda', 1e-5, 'max_iter', 10);
cfg.pilot.hp.GN    = struct('lambda', 3e-3, 'max_iter', 10);
cfg.pilot.hp.GREIT = struct('noise_figure', 0.5, 'target_size', 0.05, 'scale', 1);

%% ------------------------------------------------------------ Figuras
cfg.fig.visible = 'off';            % 'on' para verlas durante la ejecución
cfg.fig.dpi     = 200;
% Mapa de color de las imágenes de conductividad:
%   'eit'       -> fondo (dsigma = 0) azul, aumentos hacia colores cálidos
%                  (cian, verde, amarillo, rojo) y disminuciones hacia azul oscuro
%   'diverging' -> azul-blanco-rojo (blanco = sin cambio)
cfg.fig.colormap = 'eit';
end
