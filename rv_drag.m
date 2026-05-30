function [a_decay_km, a_final] = rv_drag(Nrev)

%% ----- Constants -----
mu = 398600;         % km^3/s^2
Re = 6378.1363;      % km 
omega_E = 7.2921159e-5; % rad/s (Earth rotation rate)
A = 1.3;        %m^2
m = 170;        %kg
Cd = 2;

%% ----- INPUT: initial state -----
r0 = [4771.810686400000122, -3024.529532699999891, 3804.492539700000179];   % km
v0 = [2.888529122200000, 6.849724788200001,  1.814146876000000];             % km/s

fprintf('Input r0 = [%.6f %.6f %.6f] km\n', r0);
fprintf('Input v0 = [%.6f %.6f %.6f] km/s\n\n', v0);

%% ----- Compute classical orbital elements from r0,v0 -----
r = r0(:); v = v0(:);
r_mag = norm(r); v_mag = norm(v);
h = cross(r,v); h_mag = norm(h);
e_vec = (1/mu)*((v_mag^2 - mu/r_mag)*r - dot(r,v)*v);
e = norm(e_vec);
energy = v_mag^2/2 - mu/r_mag;
a = -mu/(2*energy);
p = h_mag^2/mu;
i = acos(h(3)/h_mag);

% RAAN
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

% True anomaly
if e > 1e-12
    nu0 = acos(dot(e_vec,r)/(e*r_mag));
    if dot(r,v) < 0, nu0 = 2*pi - nu0; end
else
    nu0 = 0;
end

% Mean motion and period
n0 = sqrt(mu / a^3);
T0 = 2*pi / n0;

fprintf('Elements from r0,v0: a=%.6f km, e=%.6f, i=%.6f deg\n', a, e, rad2deg(i));
fprintf(' RAAN=%.6f deg, omega=%.6f deg, nu0=%.6f deg\n\n', rad2deg(RAAN), rad2deg(omega), rad2deg(nu0));
fprintf('T0 = %.6f s \n\n', T0); 

%% ----- Orbit propagation (Keplerian) -----
Nrev = 2;
if isempty(Nrev), Nrev = 1; end
dnu = 5; % degree step

nu_table = 0:dnu:360*Nrev;
N = numel(nu_table);

r_set = zeros(N,3);
v_set = zeros(N,3);
t_set = zeros(N,1);

% Rotation matrix (constant for Keplerian)
R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
R3_w = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
Q = R3_W * R1_i * R3_w;

for k = 1:N
    % Compute true anomaly in radians
    nu_rad = mod(nu0 + deg2rad(nu_table(k)), 2*pi);
    
    % Eccentric anomaly
    denom = 1 + e*cos(nu_rad);
    sinE = sqrt(1 - e^2)*sin(nu_rad)/denom;
    cosE = (e + cos(nu_rad))/denom;
    E = atan2(sinE, cosE);
    if E < 0, E = E + 2*pi; end

    % Mean anomaly and time
    M = E - e*sin(E);
    revolutions = floor(nu_table(k)/360);
    t = (M - (nu0 - e*sin(acos((e + cos(nu0))/(1 + e*cos(nu0)))))) / n0 + revolutions*T0;

    % Perifocal coordinates
    r_pf = (p/(1 + e*cos(nu_rad)))*[cos(nu_rad); sin(nu_rad); 0];
    v_pf = sqrt(mu/p)*[-sin(nu_rad); e + cos(nu_rad); 0];

    % Transform to ECI
    r_set(k,:) = (Q * r_pf)';
    v_set(k,:) = (Q * v_pf)';
    t_set(k) = t;

    % Force exact repeat at each revolution
    if mod(nu_table(k),360) == 0
        r_set(k,:) = r0;
        v_set(k,:) = v0;
        t_set(k) = revolutions*T0;
    end
end

%% ----- Convert ECI → ECEF → Geodetic (lat, lon, alt) -----
lat_deg = zeros(N,1);
lon_deg = zeros(N,1);
alt_km = zeros(N,1);

for k = 1:N
    % Earth rotation angle
    theta_g = omega_E * t_set(k);

    % ECI → ECEF rotation
    R3 = [cos(theta_g) sin(theta_g) 0; -sin(theta_g) cos(theta_g) 0; 0 0 1];
    r_ecef = R3 * r_set(k,:)';

    x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);
    r_xy = sqrt(x^2 + y^2);

    % Latitude 
    lat_rad = atan2(z, r_xy);
    lon_rad = atan2(y, x);

    lat_deg(k) = rad2deg(lat_rad);
    lon_deg(k) = rad2deg(lon_rad);
    alt_km(k) = norm(r_ecef) - Re;
end

%% ----- Atmospheric Model + Drag + Semi-major Axis Decay -----
year = 2025; day_of_year = 329; UTseconds = 0;
f107A = 100; f107 = 100; ap = [4 0 0 0 0 0 0];

temp_array = zeros(1,N);
density_array = zeros(1,N);
da_dt_array = zeros(1,N);

% initialize a(t) in meters
a_decay = zeros(1,N);
a_decay(1) = a * 1000;      % convert km → m

for k = 1:N

    % 1. Atmospheric model 
    alt_m = alt_km(k) * 1000;    
    [T, rho] = atmosnrlmsise00(alt_m, lat_deg(k), lon_deg(k), ...
                               year, day_of_year, UTseconds, ...
                               f107A, f107, ap);

    temp_array(k)    = T(2);
    density_array(k) = rho(6);

    % 2. Relative velocity
    theta_g = omega_E * t_set(k);
    R3 = [cos(theta_g) sin(theta_g) 0;
         -sin(theta_g) cos(theta_g) 0; 
          0            0           1];

    v_eci_m = v_set(k,:)' * 1000;
    r_eci_m = r_set(k,:)' * 1000;

    v_ecef  = R3 * v_eci_m;
    r_ecef  = R3 * r_eci_m;

    v_atm = cross([0;0;omega_E], r_ecef);
    v_rel = norm(v_ecef - v_atm);

    % 3. Mean motion using CURRENT a(t)
    a_curr_km = a_decay(k) / 1000;     % convert m → km
    n = sqrt(mu / a_curr_km^3);        % rad/s

    % 4. True anomaly scalar at point k
    nu_rad = mod(nu0 + deg2rad(nu_table(k)), 2*pi);

    % 5. Drag decay rate da/dt using updated a(t)
    sqrt_term  = sqrt(max(0, 1 + e^2 + 2*e*cos(nu_rad)));
    denom_term = n * sqrt(max(0, 1 - e^2));

    da_dt = - density_array(k) * (Cd * A / m) * ...
             (v_rel^2) * (sqrt_term / denom_term);   % m/s

    da_dt_array(k) = da_dt;

    % 6. Integrate to next step
    if k < N
        dt = t_set(k+1) - t_set(k);
        a_decay(k+1) = a_decay(k) + da_dt * dt;
    end
end

a_decay_km = a_decay / 1000;   % convert to km

%% ----- Display result -----
idx = round(N);
fprintf('\nSample point at t = %.15f s:\n', t_set(idx));
fprintf('  Latitude  = %.15f°\n', lat_deg(idx));
fprintf('  Longitude = %.15f°\n', lon_deg(idx));
fprintf('  Altitude  = %.15f km\n', alt_km(idx));
fprintf('  Temperature = %.15f K\n', temp_array(idx));
fprintf('  Density     = %.15e kg/m^3\n', density_array(idx));
fprintf('  da/dt       = %.15e m/s\n', da_dt_array(idx));
fprintf('  a (initial) = %.15f km\n', a);
fprintf('  a (decayed) = %.15f km\n', a_decay_km(idx));

%% ----- Plot Semi-major Axis Decay (with reference and end difference) -----

days = 1;

% Reference semi-major axis (ideal)
a_ref = 7119.041210;   % km

% create smooth time grid (seconds)
t_query_a = linspace(0, days*86400, 2000);

% interpolate a(t) to user-defined grid
a_interp = interp1(t_set, a_decay_km, t_query_a, 'linear', 'extrap');
a_final  = a_interp(end);

figure; hold on; grid on; box on;

% Decayed semi-major axis
plot(t_query_a/86400, a_interp, 'b-', 'LineWidth', 2);

% Reference (constant) semi-major axis
plot(t_query_a/86400, a_ref*ones(size(t_query_a)), 'k--', 'LineWidth', 1.5);

% End-day marker
t_end = t_query_a(end)/86400;
plot(t_end, a_final, 'bo', 'MarkerFaceColor','b', 'MarkerSize',7);

% End-day difference
da_end = a_final - a_ref;

% Annotate ONLY end difference
text(t_end, a_final, ...
    sprintf('  \\Delta a = %.4f km', da_end), ...
    'FontSize',10, 'FontWeight','bold', 'Color','b');

xlabel('Time (days)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Semi-major Axis a(t) (km)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Semi-major Axis Decay due to Atmospheric Drag (%.0f day)', days), ...
      'FontSize', 13, 'FontWeight', 'bold');

legend({'Decayed semi-major axis a(t)', ...
        'Ideal orbit semi-major axis a_0'}, ...
        'Location','best');