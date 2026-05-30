%% ----- Constants -----
mu = 398600;         % km^3/s^2
Re = 6378.1363;      % km 

%% ----- INPUT: initial state (ECI) -----
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

% True anomaly at epoch
if e > 1e-12
    nu0 = acos(dot(e_vec,r)/(e*r_mag));
    if dot(r,v) < 0, nu0 = 2*pi - nu0; end
else
    nu0 = 0;
end

% Mean motion and period
n0 = sqrt(mu / a^3);
T0 = 2*pi / n0;

fprintf('a = %.6f km, e = %.6f, i = %.6f deg\n', a, e, rad2deg(i));
fprintf('RAAN = %.6f deg, omega = %.6f deg, nu0 = %.6f deg\n', rad2deg(RAAN), rad2deg(omega), rad2deg(nu0));
fprintf('Period = %.6f s (%.6f hr)\n\n', T0, T0/3600);

%% ----- Rotation matrix -----
R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
R3_w = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
Q = R3_W * R1_i * R3_w;

%% ----- Orbit propagation -----
Nrev = input('Enter number of revolutions: ');
if isempty(Nrev), Nrev = 1; end

dnu = 5;                      % only for desired resolution
nu_table = 0:dnu:360*Nrev;    % just used for length / labeling
N = numel(nu_table);

% uniform time samples over Nrev periods
t_set = linspace(0, Nrev*T0, N).';   % s

r_set = zeros(N,3);
v_set = zeros(N,3);
nu_set = zeros(N,1);

% Initial mean anomaly
E0 = 2 * atan2( sqrt(1-e) * sin(nu0/2), sqrt(1+e) * cos(nu0/2) );
M0 = E0 - e*sin(E0);

for k = 1:N
    t = t_set(k);

    % Mean anomaly at time t
    M = M0 + n0 * t;
    M = mod(M, 2*pi);

    % Solve Kepler's equation for E(t) by fixed-point iteration
    E = M;
    for iter = 1:20
        E = M + e*sin(E);
    end

    % True anomaly at time t
    nu = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    nu_set(k) = nu;

    % Perifocal position and velocity
    r_pf = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

    % Transform to ECI
    r_set(k,:) = (Q * r_pf).';
    v_set(k,:) = (Q * v_pf).';
end

%% ----- Output table -----
fprintf('\n%6s %15s %15s %15s %15s %15s %15s %15s\n', ...
    'ν(deg)','t(s)','r_x(km)','r_y(km)','r_z(km)', ...
    'v_x(km/s)','v_y(km/s)','v_z(km/s)');
fprintf('%s\n', repmat('-',1,140));

for k = 1:N
    fprintf('%6.0f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n', ...
        nu_table(k), t_set(k), ...
        r_set(k,1), r_set(k,2), r_set(k,3), ...
        v_set(k,1), v_set(k,2), v_set(k,3));
end

%% ----- Plot orbit -----
figure; hold on; grid on; axis equal;
plot3(r_set(:,1), r_set(:,2), r_set(:,3), 'b-', 'LineWidth', 1.5);

% Earth (scaled)
[xs, ys, zs] = sphere(60);
surf(0.4*Re*xs, 0.4*Re*ys, 0.4*Re*zs, ...
     'FaceColor',[0.2 0.5 1], 'EdgeColor','none', 'FaceAlpha',0.3);

% Markers
rmag = vecnorm(r_set,2,2);
[~, idx_perigee] = min(rmag);
[~, idx_apogee] = max(rmag);

plot3(r_set(1,1),          r_set(1,2),          r_set(1,3),          'go', 'MarkerFaceColor','g'); % initial
plot3(r_set(idx_perigee,1),r_set(idx_perigee,2),r_set(idx_perigee,3),'ro', 'MarkerFaceColor','r'); % perigee
plot3(r_set(idx_apogee,1), r_set(idx_apogee,2), r_set(idx_apogee,3), 'mo', 'MarkerFaceColor','m'); % apogee

days_query = input('Enter time in DAYS: ');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('Keplerian Orbit - %d days', days_query));
legend('Orbit','Earth','Initial','Perigee','Apogee','Location','bestoutside');
view(35,25);

%% ----- Query state at user-defined time (DAYS) -----
t_query = days_query * 86400;

%% ----- Compute Mean Anomaly at t_query -----
M_t = M0 + n0 * t_query;
M_t = mod(M_t, 2*pi);

%% ----- Solve Kepler's Equation for E(t_query) -----
E = M_t;
for iter = 1:20
    E = M_t + e*sin(E);
end

%% ----- True anomaly at time t_query -----
nu_t = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

%% ----- Perifocal position and velocity at t_query -----
r_pf = (p / (1 + e*cos(nu_t))) * [cos(nu_t); sin(nu_t); 0];
v_pf = sqrt(mu/p) * [-sin(nu_t); e + cos(nu_t); 0];

%% ----- Convert to ECI -----
r_tq = Q * r_pf;
v_tq = Q * v_pf;

%% ----- Print results -----
fprintf('\nState Vector at t = %.0f days:\n', days_query);
fprintf('r(t) = [%.6f  %.6f  %.6f] km\n', r_tq(1), r_tq(2), r_tq(3));
fprintf('v(t) = [%.6f  %.6f  %.6f] km/s\n\n', v_tq(1), v_tq(2), v_tq(3));

%% ----- Mark the point on the orbit plot -----
figure(1); hold on;
plot3(r_tq(1), r_tq(2), r_tq(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor','k');
text(r_tq(1), r_tq(2), r_tq(3), sprintf('  t=%.0f days', days_query), ...
     'Color','k', 'FontSize', 9);
title(sprintf('Keplerian Orbit - %d days', days_query));
legend('Orbit','Earth','Initial','Perigee','Apogee','Input day','Location','bestoutside');