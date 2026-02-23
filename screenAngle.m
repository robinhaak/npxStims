clc; clear; close all;

% =========================================================
% MODE: 'forward'  -> give r1, r2, S, theta_near => get eccentricities + plot
%       'inverse'  -> give ecc_near, ecc_far, S   => get r1, r2, theta_near + plot
% =========================================================
mode = 'forward'; %'inverse';

eye           = [0, 0];
mouse_forward = [0, 1];

% --- FORWARD inputs ---
r1        = 36; %24;       % cm
r2        = 38; %53;       % cm
S         = 51;       % cm
theta_near = 3*pi/4;  % rad, direction to near edge from x-axis

% --- INVERSE inputs ---
ecc_near_target = -64;   % deg, right of forward (near edge) %from Jorrits paper
ecc_far_target  = +25;   % deg, left  of forward (far edge)
S_inv           = 51;    % cm, screen physical length
r1_inv          = 24;    % cm, fix distance to near edge, r2 will be solved

% =========================================================

fwd_angle = atan2(mouse_forward(2), mouse_forward(1));

if strcmp(mode, 'inverse')
    S  = S_inv;
    r1 = r1_inv;

    % Direction angles from x-axis
    angle_near = fwd_angle - deg2rad(ecc_near_target);
    angle_far  = fwd_angle - deg2rad(ecc_far_target);

    u_near = [cos(angle_near), sin(angle_near)];
    u_far  = [cos(angle_far),  sin(angle_far)];

    % Law of cosines: S^2 = r1^2 + r2^2 - 2*r1*r2*cos(theta)
    cos_theta = dot(u_near, u_far);
    % Solve quadratic for r2
    A =  1;
    B = -2*r1*cos_theta;
    C =  r1^2 - S^2;
    disc = B^2 - 4*A*C;
    if disc < 0
        error('No real solution for r2 ??? check inputs');
    end
    r2_candidates = [(-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A)];
    r2 = r2_candidates(r2_candidates > 0);
    if isempty(r2)
        error('No positive r2 solution');
    end
    r2 = r2(1);  % take first positive

    near_edge  = r1 * u_near;
    far_edge   = r2 * u_far;
    theta_near = angle_near;

    fprintf('=== INVERSE SOLUTION ===\n');
    fprintf('r1 (near edge distance): %.2f cm\n', r1);
    fprintf('r2 (far edge distance):  %.2f cm\n', r2);
    fprintf('S  (screen length):      %.2f cm\n', S);
    fprintf('theta_near:              %.4f rad  (%.1f deg from x-axis)\n', theta_near, rad2deg(theta_near));
    fprintf('Ecc near: %.2f deg  |  Ecc far: %.2f deg\n', ecc_near_target, ecc_far_target);

else  % forward
    near_edge = r1 * [cos(theta_near), sin(theta_near)];

    d = norm(near_edge - eye);
    if d > r2+S || d < abs(r2-S)
        error('No solution: check r1, r2, S');
    end
    a  = (r2^2 - S^2 + d^2) / (2*d);
    h  = sqrt(r2^2 - a^2);
    P2 = eye + a*(near_edge - eye)/d;
    perp = h * [-(near_edge(2)-eye(2))/d, (near_edge(1)-eye(1))/d];
    F1 = P2 + perp;
    F2 = P2 - perp;
    far_edge = F1;
    if F2(2) > F1(2), far_edge = F2; end
end

% --- Visual angle ---
cos_va = (r1^2 + r2^2 - S^2) / (2*r1*r2);
visual_angle_deg = rad2deg(acos(cos_va));

% --- Eccentricities ---
angle_near_c = atan2(near_edge(2)-eye(2), near_edge(1)-eye(1));
angle_far_c  = atan2(far_edge(2)-eye(2),  far_edge(1)-eye(1));
ecc_near = rad2deg(fwd_angle - angle_near_c);
ecc_far  = rad2deg(fwd_angle - angle_far_c);

fprintf('Ecc near (computed): %.2f deg\n', ecc_near);
fprintf('Ecc far  (computed): %.2f deg\n', ecc_far);
fprintf('Visual angle:        %.2f deg\n', visual_angle_deg);

% --- Rotate so forward = straight up ---
rot_angle = pi/2 - fwd_angle;
R = [cos(rot_angle), -sin(rot_angle);
     sin(rot_angle),  cos(rot_angle)];

near_r = (R * near_edge')';
far_r  = (R * far_edge')';
eye_r  = (R * eye')';

% --- Derived geometry in rotated frame ---
screen_dir = far_r - near_r;

% Perpendicular foot from eye to screen
t_perp    = dot(eye_r - near_r, screen_dir) / dot(screen_dir, screen_dir);
perp_foot = near_r + t_perp * screen_dir;
perp_dist = norm(perp_foot - eye_r);

% Forward-screen intersection
fwd_intersect = [];
if abs(screen_dir(1)) > 1e-9
    t_int = (eye_r(1) - near_r(1)) / screen_dir(1);
    if t_int >= 0 && t_int <= 1
        candidate = near_r + t_int * screen_dir;
        if candidate(2) > eye_r(2)
            fwd_intersect = candidate;
        end
    end
end

% Rotated angles for arcs
angle_near_r = atan2(near_r(2)-eye_r(2), near_r(1)-eye_r(1));
angle_far_r  = atan2(far_r(2)-eye_r(2),  far_r(1)-eye_r(1));
fwd_angle_r  = pi/2;

% --- PLOT ---
figure; hold on; axis equal; grid on;

% Forward line
fwd_len = max(r1,r2) * 1.2;
plot([eye_r(1), eye_r(1)], [eye_r(2), eye_r(2)+fwd_len], ...
    'k-', 'LineWidth', 1.5, 'DisplayName', 'Forward');
text(eye_r(1)+1, eye_r(2)+fwd_len+1, 'Forward', 'FontSize', 11, 'Color', 'k');

% Screen
plot([near_r(1), far_r(1)], [near_r(2), far_r(2)], ...
    'b-', 'LineWidth', 3, 'DisplayName', 'Screen');

% Eye to edges
plot([eye_r(1), near_r(1)], [eye_r(2), near_r(2)], ...
    'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('r1 = %.1f cm', r1));
plot([eye_r(1), far_r(1)],  [eye_r(2), far_r(2)],  ...
    'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('r2 = %.1f cm', r2));

% Perpendicular
if t_perp >= 0 && t_perp <= 1
    plot([eye_r(1), perp_foot(1)], [eye_r(2), perp_foot(2)], ...
        'c-', 'LineWidth', 1.5, 'DisplayName', sprintf('Perp = %.1f cm', perp_dist));
    plot(perp_foot(1), perp_foot(2), 'cs', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'c', 'HandleVisibility', 'off');
    text(perp_foot(1)+1, perp_foot(2), ...
        sprintf('d_{\\perp} = %.1f cm', perp_dist), 'FontSize', 10, 'Color', 'c');
end

% Forward-screen intersection
if ~isempty(fwd_intersect)
    plot(fwd_intersect(1), fwd_intersect(2), 'kx', ...
        'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Fwd \cap Screen');
    text(fwd_intersect(1)+1, fwd_intersect(2), ...
        'Fwd \cap Screen', 'FontSize', 10, 'Color', 'k');
end

% Eye
plot(eye_r(1), eye_r(2), 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Eye');

% Arcs
arc_r = min(r1,r2) * 0.25;

arc1 = linspace(fwd_angle_r, angle_near_r, 100);
plot(eye_r(1) + arc_r*cos(arc1), eye_r(2) + arc_r*sin(arc1), ...
    'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Ecc near = %.1f??', ecc_near));
text(eye_r(1) + arc_r*cos(mean(arc1)) - 2, ...
     eye_r(2) + arc_r*sin(mean(arc1)), ...
     sprintf('%.1f??', ecc_near), 'FontSize', 11, 'Color', 'r');

arc2 = linspace(fwd_angle_r, angle_far_r, 100);
plot(eye_r(1) + arc_r*1.5*cos(arc2), eye_r(2) + arc_r*1.5*sin(arc2), ...
    'g-', 'LineWidth', 1.5, 'DisplayName', sprintf('Ecc far = %.1f??', ecc_far));
text(eye_r(1) + arc_r*1.5*cos(mean(arc2)) - 2, ...
     eye_r(2) + arc_r*1.5*sin(mean(arc2)), ...
     sprintf('%.1f??', ecc_far), 'FontSize', 11, 'Color', 'g');

arc3 = linspace(angle_near_r, angle_far_r, 100);
plot(eye_r(1) + arc_r*2.1*cos(arc3), eye_r(2) + arc_r*2.1*sin(arc3), ...
    'm-', 'LineWidth', 2, 'DisplayName', sprintf('Visual angle = %.1f??', visual_angle_deg));
text(eye_r(1) + arc_r*2.1*cos(mean(arc3)) - 2, ...
     eye_r(2) + arc_r*2.1*sin(mean(arc3)), ...
     sprintf('%.1f??', visual_angle_deg), 'FontSize', 13, 'Color', 'm');

% Labels
text(near_r(1)-2, near_r(2)-2, sprintf('Near edge\n(ecc=%.1f??)', ecc_near), 'FontSize', 10);
text(far_r(1)-2,  far_r(2)+2,  sprintf('Far edge\n(ecc=%.1f??)', ecc_far),  'FontSize', 10);

title(sprintf('Screen placement | S=%.0fcm | ecc: %.0f?? to %.0f??', S, ecc_near, ecc_far));
xlabel('cm'); ylabel('cm');
legend('Location', 'best');