clc; clear; close all;

%% ----- Constants -----
mu = 398600;         % km^3/s^2
Re = 6378.1363;      % km
J2 = 0.00108263;

%% ----- INPUT: initial state (ECI) -----
r0 = [4771.810686400000122, -3024.529532699999891, 3804.492539700000179];   % km
v0 = [2.888529122200000, 6.849724788200001,  1.814146876000000];             % km/s

fprintf('Input r0 = [%.6f %.6f %.6f] km\n', r0);
fprintf('Input v0 = [%.6f %.6f %.6f] km/s\n\n', v0);

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

% true anomaly at epoch (nu0)
if e > 1e-12
    nu0 = acos(dot(e_vec,r)/(e*r_mag));
    if dot(r,v) < 0, nu0 = 2*pi - nu0; end
else
    nu0 = 0;
end

% mean motions / periods
n0 = sqrt(mu / a^3);           % Keplerian mean motion (rad/s)
T0 = 2*pi / n0;

% J2-corrected mean motion
n_J2 = n0 * (1 + 1.5*J2*(Re^2/p^2)*sqrt(1 - e^2)*(1 - 1.5*sin(i)^2));
T_J2 = 2*pi / n_J2;

fprintf('Elements from r0,v0: a=%.6f km, e=%.6f, i=%.6f deg\n', a, e, rad2deg(i));
fprintf(' RAAN=%.6f deg, omega=%.6f deg, nu0=%.6f deg\n\n', rad2deg(RAAN), rad2deg(omega), rad2deg(nu0));
fprintf('T0 = %.6f s, T_J2 = %.6f s \n\n', T0, T_J2);

%% ----- J2 secular rates (rad/s) -----
Omega_dot = -1.5 * J2 * (Re^2 / p^2) * n_J2 * cos(i);                       % dRAAN/dt
omega_dot = 0.75 * J2 * (Re^2 / p^2) * n_J2 * (5*cos(i)^2 - 1);             % domega/dt

fprintf('J2 rates: dRAAN/dt = %.6e rad/s (%.6f deg/day)\n', Omega_dot, rad2deg(Omega_dot)*86400);
fprintf('         domega/dt = %.6e rad/s (%.6f deg/day)\n\n', omega_dot, rad2deg(omega_dot)*86400);

%% Precompute epoch E0 and M0 robustly from nu0
E0 = 2 * atan2( sqrt(max(0,1-e))*sin(nu0/2), sqrt(1+e)*cos(nu0/2) );
E0 = mod(E0, 2*pi);
M0 = E0 - e*sin(E0);
M0 = mod(M0, 2*pi);

%% ----- User input: how many revolutions to output ----- 
Nrev = input('Enter number of revolutions: ');
if isempty(Nrev), Nrev = 1; end

dnu = 5;  % deg step

%% 5° table steps
nu_table = 0 : dnu : 360*Nrev;
N = numel(nu_table);

% reallocate
r_set = zeros(N,3);
v_set = zeros(N,3);
t_set = zeros(N,1);
nu_out_deg = zeros(N,1);    % corrected true anomaly (deg)
M_out_deg = zeros(N,1);     % mean anomaly (deg)
RAAN_out_deg = zeros(N,1);
omega_out_deg = zeros(N,1);

for k = 1:N
    % table anomaly (0..360*Nrev)
    nu_tbl_deg = nu_table(k);

    % Absolute orbital anomaly = initial nu0 + table anomaly
    nu_total_deg = rad2deg(nu0) + nu_tbl_deg;

    % Determine completed revolutions from the absolute orbital anomaly
    rev = floor(nu_total_deg / 360);

    % Local true anomaly within the current orbit [0,360)
    nu_local_deg = mod(nu_total_deg, 360);
    nu_rad = deg2rad(nu_local_deg);

    % --- conversion of true anomaly -> eccentric anomaly
    denom = 1 + e * cos(nu_rad);
    sinE = sqrt(max(0,1 - e^2)) * sin(nu_rad) / denom;
    cosE = (e + cos(nu_rad)) / denom;
    E_from_nu = atan2(sinE, cosE);
    E_from_nu = mod(E_from_nu, 2*pi);

    t_total = (nu_tbl_deg/360) * T_J2;

    % J2 secular drift applied at t_total (same as original)
    RAAN_t = RAAN + Omega_dot * t_total;
    omega_t = omega + omega_dot * t_total;

    % Perifocal r,v for local nu_rad (same as original)
    r_pf = (p / (1 + e*cos(nu_rad))) * [cos(nu_rad); sin(nu_rad); 0];
    v_pf = sqrt(mu/p) * [-sin(nu_rad); e + cos(nu_rad); 0];

    % Perifocal -> ECI (use updated RAAN_t & omega_t) (same as original)
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
        t_total = 0;                      % define t=0 as the epoch of r0,v0
        t_since_peri = 0;
        rev = 0;
        RAAN_t = RAAN;                    % make RAAN_t equal to RAAN at epoch
        omega_t = omega;                  % make omega_t equal to omega at epoch
    end 
    
    % ----------------- NEW: compute M & nu consistent with time t_total -----
    % Propagate mean anomaly from epoch using J2-corrected mean motion
    M_t = mod(M0 + n_J2 * t_total, 2*pi);

    % Solve Kepler's equation (M_t -> E_t) with Newton
    E = M_t;
    for iter = 1:100
        E_next = E - (E - e*sin(E) - M_t) / (1 - e*cos(E));
        if abs(E_next - E) < 1e-12, break; end
        E = E_next;
    end
    E_t = mod(E,2*pi);

    % True anomaly from E_t
    nu_t = 2 * atan2( sqrt(1+e) * sin(E_t/2), sqrt(1-e) * cos(E_t/2) );
    nu_t = mod(nu_t, 2*pi);

    % ----------------- Store results -----------------
    r_set(k,:) = r_ECI;
    v_set(k,:) = v_ECI;
    t_set(k) = t_total;

    % Store the corrected anomalies (deg)
    nu_out_deg(k) = rad2deg(nu_t);
    M_out_deg(k) = rad2deg(M_t);

    RAAN_out_deg(k) = rad2deg(RAAN_t);
    omega_out_deg(k) = rad2deg(omega_t);
end

%% ----- Print table -----
fprintf('\n%6s %15s %15s %15s %15s %15s %15s %15s %12s %12s\n',...
    'ν(deg)','Time(s)','r_x(km)','r_y(km)','r_z(km)','v_x(km/s)','v_y(km/s)','v_z(km/s)','M(deg)','RAAN(deg)');
fprintf('%s\n', repmat('-',1,150));

for k = 1:N
    fprintf('%6.2f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %12.6f %12.6f\n', ...
        nu_out_deg(k), t_set(k), r_set(k,1), r_set(k,2), r_set(k,3), v_set(k,1), v_set(k,2), v_set(k,3), M_out_deg(k), RAAN_out_deg(k));
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
title(sprintf('Orbit with J2 (%d revs, Δν = %d°)', Nrev, dnu));
legend('Orbit','Earth','Initial','Perigee','Apogee','Location','bestoutside');
view(35,25);

%% ----- plot secular drift of RAAN and omega vs time (user days input) -----
days = input('\nEnter number of days to see RAAN & ω drift: ');
if isempty(days), days = 2; end

% ----- Constant reference values -----
RAAN_const  = 265.771134;      % deg
omega_const = 115.120347;     % deg

t_days = linspace(0, days*86400, 3600);  % seconds
RAAN_drift  = rad2deg(RAAN + Omega_dot * t_days);
omega_drift = rad2deg(omega + omega_dot * t_days);

figure; hold on; grid on;

% J2 secular drift
plot(t_days/86400, RAAN_drift, 'b-', 'LineWidth', 1.5);
plot(t_days/86400, omega_drift, 'r-', 'LineWidth', 1.5);

% Constant reference lines
plot(t_days/86400, RAAN_const*ones(size(t_days)), 'g--', 'LineWidth', 1.2);
plot(t_days/86400, omega_const*ones(size(t_days)), 'k--', 'LineWidth', 1.2);

% ----- End-day difference -----
t_end = t_days(end)/86400;

RAAN_end  = RAAN_drift(end);
omega_end = omega_drift(end);

dRAAN_end  = RAAN_end  - RAAN_const;
domega_end = omega_end - omega_const;

% Mark end-day points
plot(t_end, RAAN_end,  'bo', 'MarkerFaceColor','b', 'MarkerSize',7);
plot(t_end, omega_end, 'ro', 'MarkerFaceColor','r', 'MarkerSize',7);

% Annotate ONLY end-day difference
text(t_end, RAAN_end, ...
    sprintf('  ΔΩ = %.4f°', dRAAN_end), ...
    'Color','b','FontSize',10,'FontWeight','bold');

text(t_end, omega_end, ...
    sprintf('  Δω = %.4f°', domega_end), ...
    'Color','r','FontSize',10,'FontWeight','bold');

xlabel('Time (days)');
ylabel('Angle (deg)');
title(sprintf('J_2 Secular Drift of RAAN and \\omega over %.0f days', days));

legend({'RAAN (Ω) – J_2 drift', ...
        'Argument of perigee (ω) – J_2 drift', ...
        'Ideal orbit RAAN', ...
        'Ideal orbit ω'}, ...
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

%% ----- Compute state vectors at user-specified time (seconds) -----
t_query = days*86400;

if ~isempty(t_query) && t_query >= 0
    % --- Interpolate orbit table for exact position ---
    if t_query <= t_set(end)
        % Interpolate position and velocity on the precomputed orbit
        r_tq = interp1(t_set, r_set, t_query, 'linear')';
        v_tq = interp1(t_set, v_set, t_query, 'linear')';
    else
        % For times beyond table, propagate using J2-corrected mean motion
        RAAN_tq = RAAN + Omega_dot * t_query;
        omega_tq = omega + omega_dot * t_query;

        % Compute mean anomaly at t_query from epoch using n_J2
        M_tq = mod(M0 + n_J2 * t_query, 2*pi);

        % Solve Kepler's equation
        E = M_tq;
        for iter = 1:200
            E_next = E - (E - e*sin(E) - M_tq) / (1 - e*cos(E));
            if abs(E_next - E) < 1e-13, break; end
            E = E_next;
        end
        E = mod(E,2*pi);

        % True anomaly
        nu_tq = 2 * atan2( sqrt(1+e) * sin(E/2), sqrt(1-e) * cos(E/2) );

        % Perifocal coordinates
        r_pf = (p / (1 + e*cos(nu_tq))) * [cos(nu_tq); sin(nu_tq); 0];
        v_pf = sqrt(mu/p) * [-sin(nu_tq); e + cos(nu_tq); 0];

        % ECI conversion
        R3_W = [cos(RAAN_tq) -sin(RAAN_tq) 0;
                sin(RAAN_tq)  cos(RAAN_tq) 0;
                0              0            1];
        R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
        R3_w = [cos(omega_tq) -sin(omega_tq) 0;
                sin(omega_tq)  cos(omega_tq) 0;
                0               0             1];
        Q = R3_W * R1_i * R3_w;

        r_tq = Q * r_pf;
        v_tq = Q * v_pf;
        M_tq_deg = rad2deg(M_tq);
    end

    % --- Compute true anomaly from interpolated position ---
    % Build rotated eccentric vector for time t_query
    RAAN_tq = RAAN + Omega_dot * t_query;
    omega_tq = omega + omega_dot * t_query;
    R3_W_t = [cos(RAAN_tq) -sin(RAAN_tq) 0; sin(RAAN_tq) cos(RAAN_tq) 0; 0 0 1];
    R1_i_t = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w_t = [cos(omega_tq) -sin(omega_tq) 0; sin(omega_tq) cos(omega_tq) 0; 0 0 1];
    Q_t = R3_W_t * R1_i_t * R3_w_t;
    e_vec_t = (Q_t * [e;0;0])';

    rmag_tq = norm(r_tq);
    if e > 1e-12
        nu_from_r = acos( max(-1, min(1, dot(e_vec_t, r_tq)/(e*rmag_tq))) );
        if dot(r_tq, v_tq) < 0
            nu_from_r = 2*pi - nu_from_r;
        end
        nu_display = rad2deg(nu_from_r);
    else
        % compute from M_tq/E
        if exist('E','var')
            nu_from_r = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
            nu_display = rad2deg(mod(nu_from_r,2*pi));
        else
            nu_display = NaN;
        end
    end

    % Compute mean anomaly for t_query consistently
    M_tq = mod(M0 + n_J2 * t_query, 2*pi);

    % --- Display ---
    fprintf('\nState at t = %.0f days (from r0,v0 epoch):\n', days);
    fprintf('  r = [%.6f, %.6f, %.6f] km\n', r_tq(1), r_tq(2), r_tq(3));
    fprintf('  v = [%.6f, %.6f, %.6f] km/s\n', v_tq(1), v_tq(2), v_tq(3));
    fprintf('  RAAN = %.6f deg\n', rad2deg(RAAN_tq));
    fprintf('  ω    = %.6f deg\n', rad2deg(omega_tq));
    fprintf('  ν    = %.6f deg\n', nu_display);
    fprintf('  a    = %.6f km\n', a);
    fprintf('  e    = %.6f\n', e);
    fprintf('  i    = %.6f deg\n', rad2deg(i));
    fprintf('  M    = %.6f deg\n', rad2deg(M_tq));

    % --- Plot on orbit ---
    figs = findall(0, 'Type', 'figure');
    if ~isempty(figs), figure(figs(end)); else 
        figure; hold on; grid on; axis equal; end
    hold on;
    plot3(r_tq(1), r_tq(2), r_tq(3), 'kp', 'MarkerFaceColor','y', 'MarkerSize',10);
    text(r_tq(1), r_tq(2), r_tq(3), sprintf('  t=%.0f days', days), ...
         'Color','k','FontSize',9,'FontWeight','bold');
    legend({'Orbit','Earth','Initial','Perigee','Apogee','Queried point'}, 'Location','bestoutside');
else
    fprintf('\nNo valid time entered — skipping state vector computation.\n');
end

%% ----- Plot new orbit passing through queried state vector -----
if exist('r_tq','var') && exist('v_tq','var')
    % Compute orbital elements from queried state
    r_vec = r_tq(:); v_vec = v_tq(:);
    r_mag_q = norm(r_vec); v_mag_q = norm(v_vec);
    h_vec = cross(r_vec, v_vec); h_mag_q = norm(h_vec);
    e_vec_q = (1/mu) * ((v_mag_q^2 - mu/r_mag_q)*r_vec - dot(r_vec,v_vec)*v_vec);
    e_q = norm(e_vec_q);
    p_q = h_mag_q^2 / mu;
    energy_q = v_mag_q^2/2 - mu/r_mag_q;
    a_q = -mu/(2*energy_q);
    i_q = acos(h_vec(3)/h_mag_q);

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
    nu_q = acos(dot(e_vec_q,r_vec)/(e_q*r_mag_q));
    if dot(r_vec,v_vec) < 0, nu_q = 2*pi - nu_q; end

    % Build true anomaly array for the new orbit
    dnu_orbit = 5;                  % step size
    nu_new_deg = 0:dnu_orbit:360;   % degrees
    N_new = length(nu_new_deg);
    r_new_set = zeros(N_new,3);

    % Shift so that the orbit starts at the queried ν
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
    figs = findall(0,'Type','figure');
    if ~isempty(figs)
        figure(figs(1)); hold on;
        plot3(r_new_set(:,1), r_new_set(:,2), r_new_set(:,3), 'r--', 'LineWidth', 1.5);
        title(sprintf('Orbit with J2 (%d revs, Δν = %d°, %d days)', Nrev, dnu, round(t_query/86400)));
        legend({'Orbit','Earth','Initial','Perigee','Apogee','Queried point','New orbit'}, 'Location','bestoutside');
    end
end