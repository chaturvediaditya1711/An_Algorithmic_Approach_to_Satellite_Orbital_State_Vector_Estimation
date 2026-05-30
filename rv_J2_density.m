clc; clear; close all;

%% ----- Constants -----
mu = 398600;         % km^3/s^2
Re = 6378.1363;      % km
J2 = 0.00108263;
omega_E = 7.2921159e-5; % rad/s (Earth rotation rate)

%% ----- INPUT: initial state -----
r0 = [4771.810686400000122, -3024.529532699999891, 3804.492539700000179];   % km
v0 = [2.888529122200000, 6.849724788200001,  1.814146876000000];             % km/s

fprintf('Input r0 = [%.6f %.6f %.6f] km\n', r0);
fprintf('Input v0 = [%.6f %.6f %.6f] km/s\n\n', v0);

Nrev = input('Enter number of revolutions: ');
%% ----- User input: how many revolutions to output ----- 
if isempty(Nrev), Nrev = 1; end

dnu = 5;  % deg step

%% 5° table steps
nu_table = 0 : dnu : 360*Nrev;
N = numel(nu_table);

%% ----- Compute classical elements from r0,v0-----
r = r0(:); v = v0(:);
r_mag = norm(r); v_mag = norm(v);
h = cross(r,v); h_mag = norm(h);
e_vec = (1/mu)*((v_mag^2 - mu/r_mag)*r - dot(r,v)*v);
e = norm(e_vec);
energy = v_mag^2/2 - mu/r_mag;
a = -mu/(2*energy);           % km
p = h_mag^2/mu;               % km
i = acos(h(3)/h_mag);         % rad

% Right Ascension of Ascending Node
n_vec = cross([0;0;1], h);
if norm(n_vec) > 0
    RAAN = acos(n_vec(1)/norm(n_vec));
    if n_vec(2) < 0, RAAN = 2*pi - RAAN; end
else
    RAAN = 0;
end

% Argument of perigee
if norm(n_vec) > 0 && e > 1e-12
    omega = acos(dot(n_vec,e_vec)/(norm(n_vec)*e));
    if e_vec(3) < 0, omega = 2*pi - omega; end
else
    omega = 0;
end

% true anomaly
if e > 1e-12
    nu0 = acos(dot(e_vec,r)/(e*r_mag));
    if dot(r,v) < 0, nu0 = 2*pi - nu0; end
else
    nu0 = 0;
end

% mean motions / periods
n0 = sqrt(mu / a^3);           % Keplerian mean motion (rad/s)
T0 = 2*pi / n0;

[a_decay_km,a_final] = rv_drag(Nrev);
idx = round(N);

p_k = a_decay_km(idx) * (1 - e^2);
% J2-corrected mean motion
n_J2 = sqrt(mu / a_decay_km(idx)^3) * (1 + 1.5*J2*(Re^2/p_k^2)*sqrt(1 - e^2)*(1 - 1.5*sin(i)^2));
T_J2 = 2*pi / n_J2;

fprintf('Elements from r0,v0: a=%.6f km, e=%.6f, i=%.6f deg\n', a, e, rad2deg(i));
fprintf(' RAAN=%.6f deg, omega=%.6f deg, nu0=%.6f deg\n\n', rad2deg(RAAN), rad2deg(omega), rad2deg(nu0));
fprintf('T0 = %.6f s, T_J2_drag = %.6f s \n\n', T0, T_J2);

%% ----- J2 secular rates (rad/s) -----
Omega_dot = -1.5 * J2 * (Re^2 / p_k^2) * n_J2 * cos(i);                       % dRAAN/dt
omega_dot = 0.75 * J2 * (Re^2 / p_k^2) * n_J2 * (5*cos(i)^2 - 1);             % domega/dt

fprintf('J2 rates: dRAAN/dt = %.15e rad/s (%.15f deg/day)\n', Omega_dot, rad2deg(Omega_dot)*86400);
fprintf('         domega/dt = %.15e rad/s (%.15f deg/day)\n\n', omega_dot, rad2deg(omega_dot)*86400);

% reallocate
r_set = zeros(N,3);
v_set = zeros(N,3);
t_set = zeros(N,1);
nu_out_deg = zeros(N,1);
RAAN_out_deg = zeros(N,1);
omega_out_deg = zeros(N,1);

% initialize mean-anomaly
E0 = 2 * atan2( sqrt(1 - e) * sin(nu0/2), sqrt(1 + e) * cos(nu0/2) );
M_prev = E0 - e*sin(E0);
t_set(1) = 0;

for k = 1:N
    % table anomaly (0..360*Nrev)
    nu_tbl_deg = nu_table(k);

    % Absolute orbital anomaly = initial nu0 + table anomaly
    nu_total_deg = rad2deg(nu0) + nu_tbl_deg;

    % Determine completed revolutions from the absolute orbital anomaly
    rev = floor(nu_total_deg / 360);           % integer number of 360° passed

    % Local true anomaly within the current orbit [0,360)
    nu_local_deg = mod(nu_total_deg, 360);
    nu_rad = deg2rad(nu_local_deg);

    % --- Stable conversion true anomaly -> eccentric anomaly using atan2 ---
    denom = 1 + e * cos(nu_rad);
    sinE = sqrt(1 - e^2) * sin(nu_rad) / denom;
    cosE = (e + cos(nu_rad)) / denom;
    E = atan2(sinE, cosE);
    E = mod(E, 2*pi);

    % Mean anomaly
    M = E - e*sin(E);
    if k == 1
        dt = 0;
    else
        dM = M - M_prev;
        if dM < 0
            dM = dM + 2*pi;
        end

        % --- Use decayed semi-major axis for this table step ---
        if numel(a_decay_km) >= k
            a_k = a_decay_km(k);
        else
            a_k = a_decay_km(end);
        end

        % local p and mean motion (Keplerian)
        p_k = a_k * (1 - e^2);
        n_k = sqrt(mu / a_k^3);

        % include first-order J2 correction in mean motion
        n_k = n_k * (1 + 1.5*J2*(Re^2/p_k^2)*sqrt(1 - e^2)*(1 - 1.5*sin(i)^2));

        % time increment from mean anomaly change
        dt = dM / n_k;
    end

    % Accumulate time
    if k == 1
        t_set(k) = 0;
    else
        t_set(k) = t_set(k-1) + dt;
    end
    M_prev = M;

    % Use local p_k and n_k to compute local J2 rates for RAAN and omega drift
    if exist('p_k','var') && exist('n_k','var')
        Omega_dot_local = -1.5 * J2 * (Re^2 / p_k^2) * n_k * cos(i);
        omega_dot_local = 0.75 * J2 * (Re^2 / p_k^2) * n_k * (5*cos(i)^2 - 1);
    else
        Omega_dot_local = Omega_dot;
        omega_dot_local = omega_dot;
    end

    t_total = t_set(k);

    % J2 secular drift applied at t_total
    RAAN_t = RAAN + Omega_dot_local * t_total;
    omega_t = omega + omega_dot_local * t_total;

    % Perifocal r,v
    if exist('p_k','var')
        r_pf = (p_k / (1 + e*cos(nu_rad))) * [cos(nu_rad); sin(nu_rad); 0];
        v_pf = sqrt(mu/p_k) * [-sin(nu_rad); e + cos(nu_rad); 0];
    else
        r_pf = (p / (1 + e*cos(nu_rad))) * [cos(nu_rad); sin(nu_rad); 0];
        v_pf = sqrt(mu/p) * [-sin(nu_rad); e + cos(nu_rad); 0];
    end

    % Perifocal -> ECI
    R3_W = [cos(RAAN_t) -sin(RAAN_t) 0;
            sin(RAAN_t)  cos(RAAN_t) 0;
            0            0           1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w = [cos(omega_t) -sin(omega_t) 0;
            sin(omega_t)  cos(omega_t) 0;
            0             0            1];
    Q = R3_W * R1_i * R3_w;

    r_ECI = (Q * r_pf)';
    v_ECI = (Q * v_pf)';

    % ----------------- CONSISTENCY FIX FOR THE FIRST ROW -----------------
    if k == 1
        % Force state vectors to exact inputs
        r_ECI = r0;
        v_ECI = v0;

        % Force the epoch/time/angles to match the computed orbital elements at epoch
        t_set(1) = 0;                      % define t=0 as the epoch of r0,v0
        M_prev = M;                        % ensure M_prev aligned with epoch
        RAAN_t = RAAN;
        omega_t = omega;
    end 
    
    % Store results
    r_set(k,:) = r_ECI;
    v_set(k,:) = v_ECI;
    nu_out_deg(k) = nu_tbl_deg;
    RAAN_out_deg(k) = rad2deg(RAAN_t);
    omega_out_deg(k) = rad2deg(omega_t);
end

%% ----- Print table -----
fprintf('\n%6s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n',...
    'ν(deg)','Time(s)','r_x(km)','r_y(km)','r_z(km)','v_x(km/s)','v_y(km/s)','v_z(km/s)','RAAN(deg)','ω(deg)');
fprintf('%s\n', repmat('-',1,200));

for k = 1:N
    fprintf('%6.0f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', ...
        nu_out_deg(k), t_set(k), r_set(k,1), r_set(k,2), r_set(k,3), v_set(k,1), v_set(k,2), v_set(k,3), RAAN_out_deg(k), omega_out_deg(k));
end

%% ----- Plot orbit -----
figure; hold on; grid on; axis equal;
plot3(r_set(:,1), r_set(:,2), r_set(:,3), 'b-', 'LineWidth', 1.5);

% Earth
scale = 0.3;
[xs, ys, zs] = sphere(80);
surf(scale*Re*xs, scale*Re*ys, scale*Re*zs, 'FaceColor',[0.2 0.5 1], 'FaceAlpha', 0.25, 'EdgeColor','none');

% mark initial, perigee, apogee
rmag = sqrt(sum(r_set.^2,2));
[~, idx_perigee] = min(rmag);
[~, idx_apogee] = max(rmag);
plot3(r_set(1,1), r_set(1,2), r_set(1,3), 'go','MarkerFaceColor','g','MarkerSize',6);
plot3(r_set(idx_perigee,1), r_set(idx_perigee,2), r_set(idx_perigee,3), 'ro','MarkerFaceColor','r','MarkerSize',8);
plot3(r_set(idx_apogee,1), r_set(idx_apogee,2), r_set(idx_apogee,3), 'mo','MarkerFaceColor','m','MarkerSize',8);

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('Orbit with J2 secular drift and Drag (%d revs, Δν = %d°)', Nrev, dnu));
legend('Orbit','Earth','Initial','Perigee','Apogee','Location','bestoutside');
view(35,25);
%% ----- plot secular drift of RAAN and omega vs time (w.r.t. constant reference) -----
days = input('\nEnter number of days for variation of orbital elements plot: ');
if isempty(days), days = 2; end

% ----- Constant reference values -----
RAAN_const  = 265.771134;      % deg
omega_const = 115.120347;     % deg

t_days = linspace(0, days*86400, 3600);  % seconds

RAAN_drift  = rad2deg(RAAN  + Omega_dot * t_days);
omega_drift = rad2deg(omega + omega_dot * t_days);

figure; hold on; grid on;

% Secular drift curves
plot(t_days/86400, RAAN_drift,  'b-', 'LineWidth', 1.5);
plot(t_days/86400, omega_drift, 'r-', 'LineWidth', 1.5);

% Ideal (constant) reference lines
plot(t_days/86400, RAAN_const  * ones(size(t_days)), 'b--', 'LineWidth', 1.2);
plot(t_days/86400, omega_const * ones(size(t_days)), 'r--', 'LineWidth', 1.2);

% ----- End-day values -----
t_end = t_days(end) / 86400;

RAAN_end  = RAAN_drift(end);
omega_end = omega_drift(end);

% Differences w.r.t. constant reference
dRAAN_end  = RAAN_end  - RAAN_const;
domega_end = omega_end - omega_const;

% Mark end-day points
plot(t_end, RAAN_end,  'bo', 'MarkerFaceColor','b', 'MarkerSize',7);
plot(t_end, omega_end, 'ro', 'MarkerFaceColor','r', 'MarkerSize',7);

% Annotate ONLY final differences
text(t_end, RAAN_end, ...
    sprintf('  \\Delta\\Omega = %.4f^\\circ', dRAAN_end), ...
    'Color','b','FontSize',10,'FontWeight','bold');

text(t_end, omega_end, ...
    sprintf('  \\Delta\\omega = %.4f^\\circ', domega_end), ...
    'Color','r','FontSize',10,'FontWeight','bold');

xlabel('Time (days)');
ylabel('Angle (deg)');
title(sprintf('J_2 Secular Drift of RAAN and \\omega at %.0f days', days));

legend({'RAAN (Omega) – J_2 drift', ...
        'ω – J_2 drift', ...
        'Ideal RAAN', ...
        'Ideal ω'}, ...
        'Location','best');

% Compute final RAAN and omega after given number of days
RAAN_final = rad2deg(RAAN + Omega_dot * t_days(end));
omega_final = rad2deg(omega + omega_dot * t_days(end));

% Display results
fprintf('\nAfter %.0f days:\n', days);
fprintf('  RAAN (Ω)  = %.6f deg\n', RAAN_final);
fprintf('  ω (arg. of perigee) = %.6f deg\n', omega_final);
fprintf('  ΔRAAN = %.6f deg/day,  Δω = %.6f deg/day\n', ...
        rad2deg(Omega_dot)*86400, rad2deg(omega_dot)*86400);

%% ----- Compute state vectors at user-specified time using decayed a_final -----
t_query = days * 86400;   % seconds

% Use fully decayed semi-major axis
a_query = a_final;          % km
p_query = a_query * (1 - e^2);
n_query = sqrt(mu / a_query^3);  % Keplerian mean motion

% Include J2 correction for mean motion
n_query = n_query * (1 + 1.5*J2*(Re^2/p_query^2)*sqrt(1 - e^2)*(1 - 1.5*sin(i)^2));

% Compute RAAN and omega at t_query
RAAN_tq = RAAN + Omega_dot * t_query;
omega_tq = omega + omega_dot * t_query;

% Compute mean anomaly at t_query
E0 = 2 * atan2( sqrt(1 - e) * sin(nu0/2), sqrt(1 + e) * cos(nu0/2) );
M0 = E0 - e*sin(E0);
M_tq = mod(M0 + n_query * t_query, 2*pi);

% Solve Kepler's equation iteratively for E
E = M_tq;
for iter = 1:100
    E_next = E - (E - e*sin(E) - M_tq) / (1 - e*cos(E));
    if abs(E_next - E) < 1e-12, break; end
    E = E_next;
end

% True anomaly
nu_tq = 2 * atan2( sqrt(1+e) * sin(E/2), sqrt(1-e) * cos(E/2) );

% Perifocal coordinates
r_pf = (p_query/(1 + e*cos(nu_tq))) * [cos(nu_tq); sin(nu_tq); 0];
v_pf = sqrt(mu/p_query) * [-sin(nu_tq); e + cos(nu_tq); 0];

% Transform to ECI
R3_W = [cos(RAAN_tq) -sin(RAAN_tq) 0;
        sin(RAAN_tq)  cos(RAAN_tq) 0;
        0              0            1];
R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
R3_w = [cos(omega_tq) -sin(omega_tq) 0;
        sin(omega_tq)  cos(omega_tq) 0;
        0              0            1];
Q = R3_W * R1_i * R3_w;

r_tq = (Q*r_pf)';
v_tq = (Q*v_pf)';

% Display
fprintf('\nQueried state vector at %.0f days using decayed a_final:\n', days);
fprintf('  r = [%.6f, %.6f, %.6f] km\n', r_tq(1), r_tq(2), r_tq(3));
fprintf('  v = [%.6f, %.6f, %.6f] km/s\n', v_tq(1), v_tq(2), v_tq(3));
fprintf('  ν    = %.6f deg\n', rad2deg(nu_tq));
fprintf('  a    = %.6f km\n', a_query);
fprintf('  e    = %.6f \n', e);
fprintf('  i    = %.6f deg\n', rad2deg(i));
fprintf('  ω    = %.6f deg\n', rad2deg(omega_tq));
fprintf('  RAAN = %.6f deg\n', rad2deg(RAAN_tq));
fprintf('  M    = %.6f \n', rad2deg(M_tq));

% --- Plot on orbit ---
figs = findall(2, 'Type', 'figure');
if ~isempty(figs), figure(figs(end)); else 
    figure; hold on; grid on; axis equal; end
hold on;
plot3(r_tq(1), r_tq(2), r_tq(3), 'kp', 'MarkerFaceColor','y', 'MarkerSize',10);
text(r_tq(1), r_tq(2), r_tq(3), sprintf(' after %.0f days', days), ...
         'Color','k','FontSize',9,'FontWeight','bold');
legend({'Orbit','Earth','Initial','Perigee','Apogee','Queried point'}, 'Location','bestoutside');

%% ----- Plot new orbit passing through queried state vector -----
if exist('r_tq','var') && exist('v_tq','var')
    % Compute orbital elements from queried state
    r_vec = r_tq(:); v_vec = v_tq(:);
    r_mag = norm(r_vec); v_mag = norm(v_vec);
    h_vec = cross(r_vec, v_vec); h_mag = norm(h_vec);
    e_vec_q = (1/mu) * ((v_mag^2 - mu/r_mag)*r_vec - dot(r_vec,v_vec)*v_vec);
    e_q = norm(e_vec_q);
    p_q = h_mag^2 / mu;
    energy_q = v_mag^2/2 - mu/r_mag;
    a_q = -mu/(2*energy_q);
    i_q = acos(h_vec(3)/h_mag);

    n_vec_q = cross([0;0;1], h_vec);
    if norm(n_vec_q) > 0
        RAAN_q = acos(n_vec_q(1)/norm(n_vec_q));
        if n_vec_q(2) < 0, RAAN_q = 2*pi - RAAN_q; end
    else
        RAAN_q = 0;
    end
    if norm(n_vec_q) > 0 && e_q > 1e-12
        omega_q = acos(dot(n_vec_q, e_vec_q)/(norm(n_vec_q)*e_q));
        if e_vec_q(3) < 0, omega_q = 2*pi - omega_q; end
    else
        omega_q = 0;
    end
    nu_q = acos(dot(e_vec_q,r_vec)/(e_q*r_mag));
    if dot(r_vec,v_vec) < 0, nu_q = 2*pi - nu_q; end

    % Build true anomaly array for the new orbit
    dnu_orbit = 5;                  % step size
    nu_new_deg = 0:dnu_orbit:360;   % degrees
    N_new = length(nu_new_deg);
    r_new_set = zeros(N_new,3);

    nu_shift_deg = rad2deg(nu_q);
    nu_new_deg_shifted = mod(nu_new_deg + nu_shift_deg, 360);

    for k = 1:N_new
        nu_rad = deg2rad(nu_new_deg_shifted(k));
        % Perifocal coordinates
        r_pf = (p_q/(1 + e_q*cos(nu_rad))) * [cos(nu_rad); sin(nu_rad); 0];

        % Transformation to ECI
        R3_W = [cos(RAAN_q) -sin(RAAN_q) 0;
                sin(RAAN_q)  cos(RAAN_q) 0;
                0            0            1];
        R1_i = [1 0 0; 0 cos(i_q) -sin(i_q); 0 sin(i_q) cos(i_q)];
        R3_w = [cos(omega_q) -sin(omega_q) 0;
                sin(omega_q)  cos(omega_q) 0;
                0             0            1];
        Q = R3_W * R1_i * R3_w;

        r_new_set(k,:) = (Q*r_pf)';
    end

    % Plot new orbit in same figure
    figs = findall(2,'Type','figure');
    if ~isempty(figs)
        figure(figs(1)); hold on;
        plot3(r_new_set(:,1), r_new_set(:,2), r_new_set(:,3), 'r--', 'LineWidth', 1.5);
        title(sprintf('Orbit with J2 secular drift and Drag (%d revs, Δν = %d°, %d days)', Nrev, dnu, t_query/86400));
        legend({'Orbit','Earth','Initial','Perigee','Apogee','Queried point','New orbit'}, 'Location','bestoutside');
    end
end