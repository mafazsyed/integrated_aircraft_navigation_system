function [x_est, b_est, Ax_f_instance, Ay_f_instance, Az_f_instance, p_f_instance, q_f_instance, r_f_instance, AoA_f_instance] = integrated_navigation(c_k, d_k, t, dt)

%% Note: Rename the function in the format of SID + your student ID as on Blackboard (e.g. if your ID is 21010000, name the function as SID21010000 and submit it as SID21010000.m)
%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   c_k: input measurements, a N x 6 matrix with rows the recordings of different time instances and the columns being [A_x_m, A_y_m, A_z_m, p_m, q_m, r_m]
%   d_k: output measurements, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_GPS_m, y_GPS_m, z_GPS_m, u_GPS_m, v_GPS_m, w_GPS_m, phi_GPS_m, theta_GPS_m, psi_GPS_m, V_TAS_m, alpha_m, beta_m]
%   t: time vector
%   dt: uniform time step size
%%%%%%%%%%%%%%%%%%%%%%%
% Output Variables
%   x_est: estimated state trajectory, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_E, y_E, z_E, u, v, w, \phi, \theta, \psi, V_{wxE}, V_{wyE}, V_{wzE}]
%   b_est: estimated state trajectory, a N x 6 matrix with rows the recordings of different time instances and the columns being [b_A_x, b_A_y, b_A_z, b_p, b_q, b_r]
%   Ax_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_x. In case of no faults detected, please return a value of 0,
%   Ay_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_y. In case of no faults detected, please return a value of 0,
%   Az_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_z. In case of no faults detected, please return a value of 0,
%   p_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of p. In case of no faults detected, please return a value of 0,
%   q_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of q. In case of no faults detected, please return a value of 0,
%   r_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of r. In case of no faults detected, please return a value of 0,
%   AoA_f_instance: a number equal to the corresponding time instances t for the first detection of faults in measurements of angle of attack sensor. In case of no faults detected, please return a value of 0. 
%%%%%%%%%%%%%%%%%%%%%%%

%% Your code here
z_k = d_k;
u_k = c_k;
s_Ax = 0.01;
s_Ay = 0.01;
s_Az = 0.01;
s_p = deg2rad(0.01);
s_q = deg2rad(0.01);
s_r = deg2rad(0.01);
s_w_b_A_x = 1;
s_w_b_A_y = 1;
s_w_b_A_z = 1;
s_w_b_p = deg2rad(1);
s_w_b_q = deg2rad(1);
s_w_b_r = deg2rad(1);
s_x_E = 5;
s_y_E = 5;
s_z_E = 10;
s_u = 0.1;
s_v = 0.1;
s_w = 0.1;
s_v_tas = 0.1;
s_phi = deg2rad(0.1);
sigma_theta = deg2rad(0.1);
s_psi = deg2rad(0.1);
sigma_alpha = deg2rad(0.1);
sigma_beta = deg2rad(0.1);
x_E = z_k(1,1);
y_E = z_k(1,2);
z_E = z_k(1,3);
u_estimate = 85;
v = 0;
w = 0;
b_A_x_estimate = 0.01;
b_A_y_estimate = 0.01;
b_A_z_estimate = 0.01;
b_p_estimate = deg2rad(0.1);
b_q_estimate = deg2rad(0.1);
b_r_estimate = deg2rad(0.1);
phi = z_k(1,7);
theta = z_k(1,8);
psi = z_k(1,9);
V_wxE = 0;
V_wyE = 0;
V_wzE = 0;
state_name = {'x_{E}', 'y_{E}', 'z_{E}', 'u', 'v', 'w', '\phi', '\theta', '\psi', 'b_{A_{x}}', 'b_{A_{y}}', 'b_{A_{z}}', 'b_{p}', 'b_{q}', 'b_{r}', 'V_{wxE}', 'V_{wyE}', 'V_{wzE}'};
output_name = {'x_{GPS}', 'y_{GPS}', 'z_{GPS}', 'u_{GPS}', 'v_{GPS}', 'w_{GPS}', '\phi_{GPS}', '\theta_{GPS}', '\psi_{GPS}', 'V_{tas}', '\alpha', '\beta'};
unit = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad', 'm/s^2', 'm/s^2', 'm/s^2', 'rad/s', 'rad/s', 'rad/s', 'm/s', 'm/s', 'm/s'};
out_unit = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad', 'm/s', 'rad', 'rad'};
stdw = [s_Ax s_Ay s_Az s_p s_q s_r s_w_b_A_x s_w_b_A_y s_w_b_A_z s_w_b_p s_w_b_q s_w_b_r];
stdv = [s_x_E s_y_E s_z_E s_u s_v s_w s_phi sigma_theta s_psi s_v_tas sigma_alpha sigma_beta];
Ex_0 = [x_E y_E z_E u_estimate v w phi theta psi b_A_x_estimate b_A_y_estimate b_A_z_estimate b_p_estimate b_q_estimate b_r_estimate V_wxE V_wyE V_wzE];
stdx_0 = [0.2 0.2 0.2 0.6 0.6 0.6 0.2 0.2 0.2 12 12 12 12 12 12 0.9 0.9 0.9];
[x_est, b_est, x_cor, innov] = runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0, true);
[~, ~, ~, innov_faulty] = runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0, false);
innov_list = cell(2, 1);
innov_list{1} = innov_faulty;
innov_list{2} = innov;
data_labels = {'A_{x}', 'A_{y}', 'A_{z}', 'p', 'q', 'r', '\alpha'};
data_faulty = [b_est innov(:,11)];
threshold_list = [50 50 20 5 1 5 0.02];
leak = std(data_faulty) - mean(data_faulty);
[g_pos_list, g_neg_list, k_alarm_pos_list, k_alarm_neg_list] = implement_cusum(data_faulty, leak, threshold_list);
time_list = get_fault_onset_times(k_alarm_pos_list, k_alarm_neg_list, data_labels);
plot_cusum(g_pos_list, g_neg_list, k_alarm_pos_list, k_alarm_neg_list, threshold_list, data_labels);
plot_estimated_states(x_cor, innov_list, t, state_name, unit, output_name, out_unit, time_list);
Ax_f_instance = time_list{1};
Ay_f_instance = time_list{2};
Az_f_instance = time_list{3};
p_f_instance = time_list{4};
q_f_instance = time_list{5};
r_f_instance = time_list{6};
AoA_f_instance = time_list{7} + 3;
end
function plot_cusum(g_pos_list, g_neg_list, k_alarm_pos_list, k_alarm_neg_list, threshold_list, data_labels)
font_size = 12;
n_rows = ceil(length(data_labels) / 3);
figure
for idx=1:length(data_labels)
    subplot(n_rows, 3, idx)
    p1=plot(g_pos_list{idx});
    hold on
    p2=plot(g_neg_list{idx});
    k_alarm_pos = k_alarm_pos_list{idx};
    k_alarm_neg = k_alarm_neg_list{idx};
    if ~isempty(k_alarm_pos) && ~isempty(k_alarm_neg)
        for i=1:length(k_alarm_pos)
            p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-100 100],'-r');
        end
        for i=1:length(k_alarm_neg)
            p4=plot([k_alarm_neg(i) k_alarm_neg(i)],[-100 100],'-b');
        end
        legend([p1,p2,p3,p4],'Positive', 'Negative','+ve Test Alarm','-ve Test Alarm','Location','best')
    elseif ~isempty(k_alarm_pos)
        for i=1:length(k_alarm_pos)
            p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-100 100],'-r');
        end
        legend([p1,p2,p3],'Positive Test', 'Negative Test', '+ve Test Alarm','Location','best')
    elseif ~isempty(k_alarm_neg)
        for i=1:length(k_alarm_neg)
            p3=plot([k_alarm_neg(i) k_alarm_neg(i)],[-100 100],'-b');
        end
        legend([p1,p2,p3],'Positive Test', 'Negative Test', '-ve Test Alarm', 'Location', 'best')
    else
        legend([p1,p2],'Positive Test', 'Negative Test', 'Location','best')
    end
    yline(-threshold_list(idx))
    yline(threshold_list(idx))
    ylim([-threshold_list(idx)-0.025 threshold_list(idx)+0.025])
    xlabel('Step', 'FontSize', font_size)
    ylabel('g_t', 'FontSize', font_size)
    title(['CUSUM: ', data_labels{idx}], 'FontSize', font_size)
end
set(gcf,'WindowState','minimized');
end
function [x_est, b_est, x_cor, innov] = runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0, cyber_attack)
Ts = dt;
N = length(t);
xhat_km1_km1 = Ex_0;
P_km1_km1 = diag(stdx_0.^2);
Q=diag(stdw.^2);
R=diag(stdv.^2);
n = length(xhat_km1_km1);
m = size(u_k, 2);
p = size(z_k, 2);
u_km1 = [zeros(1,m); u_k];
stdx_cor  = zeros(N, n);
x_cor     = zeros(N, n);
K       = cell(N, 1);
innov     = zeros(N, p);
for k=1:N
    [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
    xhat_k_km1 = x_nonlin(end,:);
    [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts);
    P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1';
    H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
    Ve = (H_k * P_k_km1 * H_k' + R);
    K_k = P_k_km1 * H_k' / Ve;
    z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]);
    xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k';
    if cyber_attack
        innov(k,:) = z_k(k,:) - z_k_km1;
        innov(k, :) = innov(k, :) / sqrtm(Ve);
        if abs(innov(k, 11)) > 3
            z_k(k, 11) = z_k_km1(11);
            innov(k, 11) = 0;
        end
    end
    I_KH = eye(n) - K_k * H_k;
    P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k';
    stdx_cor(k,:) = sqrt(diag(P_km1_km1));
    x_cor(k,:) = xhat_km1_km1;
    K{k,1} = K_k;
    innov(k,:)= z_k(k,:) - z_k_km1;
    xhat_km1_km1 = xhat_k_k;
    P_km1_km1 = P_k_k;
end
x_est = [x_cor(:,1:9) x_cor(:,16:18)];
b_est = x_cor(:,10:15);
end
function time_list = get_fault_onset_times(k_alarm_pos_list, k_alarm_neg_list, data_labels)
threshold = 10;
buffer = 1000;
for idx=1:length(data_labels)
    k_alarm_pos = k_alarm_pos_list{idx};
    k_alarm_neg = k_alarm_neg_list{idx};
    if ~isempty(k_alarm_pos) && ~isempty(k_alarm_neg)
        change_points = find(ischange(k_alarm_pos, 'linear', 'Threshold', threshold));
        if ~isempty(change_points)
            pos_test_list = [k_alarm_pos(change_points(1))];
            for i = 1:length(change_points)-1
                if k_alarm_pos(change_points(i+1)) - k_alarm_pos(change_points(i)) > buffer
                    pos_test_list = [pos_test_list; k_alarm_pos(change_points(i+1))];
                end
            end
        else
            pos_test_list = k_alarm_pos(1);
        end
        time_list{idx} = pos_test_list;
    elseif ~isempty(k_alarm_pos)
        change_points = find(ischange(k_alarm_pos, 'linear', 'Threshold', threshold));
        if ~isempty(change_points)
            pos_test_list = [k_alarm_pos(change_points(1))];
            for i = 1:length(change_points)-1
                if k_alarm_pos(change_points(i+1)) - k_alarm_pos(change_points(i)) > buffer
                    pos_test_list = [pos_test_list; k_alarm_pos(change_points(i+1))];
                end
            end
        else
            pos_test_list = k_alarm_pos(1);
        end
        time_list{idx} = pos_test_list;
    elseif ~isempty(k_alarm_neg)
        change_points = find(ischange(k_alarm_neg, 'linear', 'Threshold', threshold));
        if ~isempty(change_points)
            neg_test_list = [k_alarm_neg(change_points(1))];
            for i = 1:length(change_points)-1
                if k_alarm_neg(change_points(i+1)) - k_alarm_neg(change_points(i)) > buffer
                    neg_test_list = [neg_test_list; k_alarm_neg(change_points(i+1))];
                end
            end
        else
            neg_test_list = k_alarm_neg(1);
        end
        time_list{idx} = neg_test_list;
    else
        time_list{idx} = 0;
    end
end
end
function [g_pos_list, g_neg_list, k_alarm_pos_list, k_alarm_neg_list] = implement_cusum(x_faulty, leak, threshold_list)
faulty_skip_inds = 500;
faulty_size = size(x_faulty, 1);
faulty_columns = size(x_faulty, 2);
g_pos_list = cell(faulty_columns, 1);
g_neg_list = cell(faulty_columns, 1);
k_alarm_pos_list = cell(faulty_columns, 1);
k_alarm_neg_list = cell(faulty_columns, 1);
for idx=1:faulty_columns
    straingauge = x_faulty(faulty_skip_inds:faulty_size, idx);
    theta0 = 0;
    sigma0 = 1;
    threshold_pos = threshold_list(idx);
    threshold_neg = -threshold_list(idx);
    g_pos = 0 * straingauge;
    g_neg = 0 * straingauge;
    k_alarm_pos = [];
    k_alarm_neg = [];
    s = (straingauge - theta0) / sigma0;
    for k = 1:size(straingauge, 1) - 1
        g_pos(k+1) = g_pos(k) + s(k) - leak(idx);
        g_neg(k+1) = g_neg(k) + s(k) + leak(idx);
        if g_pos(k+1) < 0
            g_pos(k+1) = 0;
        end
        if g_pos(k+1) > threshold_pos
            k_alarm_pos = [k_alarm_pos; k+1];
            g_pos(k+1) = 0;
        end
        if g_neg(k+1) > 0
            g_neg(k+1) = 0;
        end
        if g_neg(k+1) < threshold_neg
            k_alarm_neg = [k_alarm_neg; k+1];
            g_neg(k+1) = 0;
        end
    end
    if ~isempty(k_alarm_pos)
        k_alarm_pos = k_alarm_pos + faulty_skip_inds;
    end
    if ~isempty(k_alarm_neg)
        k_alarm_neg = k_alarm_neg + faulty_skip_inds;
    end
    if isempty(k_alarm_pos)
        k_alarm_pos = 0;
    end
    if isempty(k_alarm_neg)
        k_alarm_neg = 0;
    end
    g_pos_list{idx} = g_pos;
    g_neg_list{idx} = g_neg;
    k_alarm_pos_list{idx} = k_alarm_pos;
    k_alarm_neg_list{idx} = k_alarm_neg;
end
end
function plot_estimated_states(x_cor, innov_list, t, state_name, unit, output_name, out_unit, time_list)
x_label = 'Time (s)';
y_label = 'Estimation ';
font_size = 12;
n_rows = 4;
figure
for i = 1:length(state_name)
    if i > 9 && i < 16
        continue
    elseif i >= 16
        subplot_idx = i - 6;
    else
        subplot_idx = i;
    end
    subplot(n_rows, 3, subplot_idx)
    hold on
    plot(t, x_cor(:,i), 'LineWidth', 2, 'Color', [1.0 0 0]);
    hold off
    grid on
    xlabel(x_label, 'FontSize', font_size)
    ylabel([y_label unit{i}], 'FontSize', font_size)
    title([state_name{i}], 'FontSize', font_size)
    legend('Location', 'best')
end
set(gcf,'WindowState','minimized');
n_rows = ceil((length(state_name) - 12) / 3);
figure
for i = 1:length(state_name)
    if ~(i > 9 && i < 16)
        continue
    else
        subplot_idx = i - 9;
    end
    subplot(n_rows, 3, subplot_idx)
    hold on
    plot(t, x_cor(:,i), 'LineWidth', 2, 'Color', [0 0 1.0]);
    if ~isempty(time_list)
        if time_list{i-9} ~= 0
            for j=1:length(time_list{i-9})
                if j == 1
                    xline(time_list{i-9}(j)/100, '--k', 'DisplayName', 'Fault Onset Time', 'LineWidth', 2.5, 'Color', [0 0 0]);
                else
                    xline(time_list{i-9}(j)/100, '--k', 'LineWidth', 2.5, 'Color', [0 0 0], 'HandleVisibility', 'off');
                end
            end
        end
    end
    hold off
    grid on
    xlabel(x_label, 'FontSize', font_size)
    ylabel([y_label unit{i}], 'FontSize', font_size)
    title([state_name{i}], 'FontSize', font_size)
    legend('Location', 'best')
end
set(gcf,'WindowState','minimized');
figure
i = 11;
hold on
for j = 1:2
    innov = innov_list{j};
    if j == 1
        plot(t, innov(:, 11), '-b', 'DisplayName', 'faulty \alpha', 'LineWidth', 2);
    else
        plot(t, innov(:, 11), 'DisplayName', 'corrected \alpha', 'LineWidth', 2);
    end
end
if ~isempty(time_list)
    if time_list{end} ~= 0
        xline(time_list{end}(1)/100, '--k', 'DisplayName', 'Fault Onset Time', 'LineWidth', 2.5, 'Color', [0 0 0]);
    end
end
hold off
grid on
xlabel(x_label, 'FontSize', font_size)
ylabel([y_label out_unit{i}], 'FontSize', font_size)
title(['Innovation: ' output_name{i}], 'FontSize', font_size)
legend('Location', 'best')
set(gcf,'WindowState','minimized');
end
function x_dot_vector = funcf(x_vector, c_vector, t)
x_E = x_vector(1);
y_E = x_vector(2);
z_E = x_vector(3);
u = x_vector(4);
v = x_vector(5);
w = x_vector(6);
phi = x_vector(7);
theta = x_vector(8);
psi = x_vector(9);
b_A_x = x_vector(10);
b_A_y = x_vector(11);
b_A_z = x_vector(12);
b_p = x_vector(13);
b_q = x_vector(14);
b_r = x_vector(15);
V_wxE = x_vector(16);
V_wyE = x_vector(17);
V_wzE = x_vector(18);
A_x = c_vector(1);
A_y = c_vector(2);
A_z = c_vector(3);
p = c_vector(4);
q = c_vector(5);
r = c_vector(6);
g = 9.81;
x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+V_wxE;
y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+V_wyE;
z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+V_wzE;
u_dot =(A_x-b_A_x)-g*sin(theta)+(r-b_r)*v-(q-b_q)*w;
v_dot =(A_y-b_A_y)+g*cos(theta)*sin(phi)+(p-b_p)*w-(r-b_r)*u;
w_dot =(A_z-b_A_z)+g*cos(theta)*cos(phi)+(q-b_q)*u-(p-b_p)*v;
phi_dot =(p-b_p)+(q-b_q)*sin(phi)*tan(theta)+(r-b_r)*cos(phi)*tan(theta);
theta_dot =(q-b_q)*cos(phi)-(r-b_r)*sin(phi);
psi_dot =(q-b_q)*sin(phi)/cos(theta)+(r-b_r)*cos(phi)/cos(theta);
b_A_x_dot =0;
b_A_y_dot =0;
b_A_z_dot =0;
b_p_dot =0;
b_q_dot =0;
b_r_dot =0;
V_wxE_dot =0;
V_wyE_dot =0;
V_wzE_dot =0;
x_dot_vector = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot b_A_x_dot b_A_y_dot b_A_z_dot b_p_dot b_q_dot b_r_dot V_wxE_dot V_wyE_dot V_wzE_dot]';
end
function d_vector = funch(x_vector, c_vector, t)
x_E = x_vector(1);
y_E = x_vector(2);
z_E = x_vector(3);
u = x_vector(4);
v = x_vector(5);
w = x_vector(6);
phi = x_vector(7);
theta = x_vector(8);
psi = x_vector(9);
b_A_x = x_vector(10);
b_A_y = x_vector(11);
b_A_z = x_vector(12);
b_p = x_vector(13);
b_q = x_vector(14);
b_r = x_vector(15);
V_wxE = x_vector(16);
V_wyE = x_vector(17);
V_wzE = x_vector(18);
x_GPS = x_E;
y_GPS = y_E;
z_GPS = z_E;
u_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+V_wxE;
v_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+V_wyE;
w_GPS = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+V_wzE;
phi_GPS = phi;
theta_GPS = theta;
psi_GPS = psi;
V_tas = sqrt(u^2+v^2+w^2);
alpha = atan(w/u);
beta = atan(v/sqrt(u^2+w^2));
d_vector = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS phi_GPS theta_GPS psi_GPS V_tas alpha beta];
end
function [Phi, Gamma] = funcLinDisDyn(x_vector, c_m_vector, Ts)
x_E = x_vector(1);
y_E = x_vector(2);
z_E = x_vector(3);
u = x_vector(4);
v = x_vector(5);
w = x_vector(6);
phi = x_vector(7);
theta = x_vector(8);
psi = x_vector(9);
b_A_x = x_vector(10);
b_A_y = x_vector(11);
b_A_z = x_vector(12);
b_p = x_vector(13);
b_q = x_vector(14);
b_r = x_vector(15);
V_wxE = x_vector(16);
V_wyE = x_vector(17);
V_wzE = x_vector(18);
A_x_m = c_m_vector(1);
A_y_m = c_m_vector(2);
A_z_m = c_m_vector(3);
p_m = c_m_vector(4);
q_m = c_m_vector(5);
r_m = c_m_vector(6);
A = [0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                                   -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 1, 0, 0;
    0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                                   -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 0, 1, 0;
    0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                                             - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 1;
    0, 0, 0,                   0,                                        r_m - b_r,                                        b_q - q_m,                                                                                  0,                                                                             -(981*cos(theta))/100,                                                                                                     0, -1,  0,  0,  0,                    w,                   -v, 0, 0, 0;
    0, 0, 0,           b_r - r_m,                                                0,                                        p_m - b_p,                                                      (981*cos(phi)*cos(theta))/100,                                                                    -(981*sin(phi)*sin(theta))/100,                                                                                                     0,  0, -1,  0, -w,                    0,                    u, 0, 0, 0;
    0, 0, 0,           q_m - b_q,                                        b_p - p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                                    -(981*cos(phi)*sin(theta))/100,                                                                                                     0,  0,  0, -1,  v,                   -u,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                  sin(phi)*tan(theta)*(b_r - r_m) - cos(phi)*tan(theta)*(b_q - q_m),               - cos(phi)*(b_r - r_m)*(tan(theta)^2 + 1) - sin(phi)*(b_q - q_m)*(tan(theta)^2 + 1),                                                                                                     0,  0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta), 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                        cos(phi)*(b_r - r_m) + sin(phi)*(b_q - q_m),                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,            -cos(phi),             sin(phi), 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,              (sin(phi)*(b_r - r_m))/cos(theta) - (cos(phi)*(b_q - q_m))/cos(theta), - (cos(phi)*sin(theta)*(b_r - r_m))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_q - q_m))/cos(theta)^2,                                                                                                     0,  0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta), 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
    0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0];
G = [0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    -1,  0,  0,  0,                    w,                   -v,    1,    0,    0,  0,                  -w,                   v;
    0, -1,  0, -w,                    0,                    u,    0,    1,    0,  w,                   0,                  -u;
    0,  0, -1,  v,                   -u,                    0,    0,    0,    1, -v,                   u,                   0;
    0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),    0,    0,    0,  1, sin(phi)*tan(theta), cos(phi)*tan(theta);
    0,  0,  0,  0,            -cos(phi),             sin(phi),    0,    0,    0,  0,            cos(phi),           -sin(phi);
    0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),    0,    0,    0,  0, sin(phi)/cos(theta), cos(phi)/cos(theta);
    0,  0,  0,  0,                    0,                    0, 1/Ts,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0, 1/Ts,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0, 1/Ts,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0];
G = [0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
    -1,  0,  0,  0,                    w,                   -v,    1,    0,    0,    0,                  -w,                   v;
    0, -1,  0, -w,                    0,                    u,    0,    1,    0,    w,                   0,                  -u;
    0,  0, -1,  v,                   -u,                    0,    0,    0,    1,   -v,                   u,                   0;
    0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),    0,    0,    0,    1, sin(phi)*tan(theta), cos(phi)*tan(theta);
    0,  0,  0,  0,            -cos(phi),             sin(phi),    0,    0,    0,    0,            cos(phi),           -sin(phi);
    0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),    0,    0,    0,    0, sin(phi)/cos(theta), cos(phi)/cos(theta);
    0,  0,  0,  0,                    0,                    0, 1/Ts,    0,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0, 1/Ts,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0, 1/Ts,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0, 1/Ts,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                1/Ts,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                1/Ts;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
    0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0];
[Phi, Gamma] = c2d(A, G, Ts);
end
function H = funcLinDisObs(x_vector, c_m_vector, t)
x_E = x_vector(1);
y_E = x_vector(2);
z_E = x_vector(3);
u = x_vector(4);
v = x_vector(5);
w = x_vector(6);
phi = x_vector(7);
theta = x_vector(8);
psi = x_vector(9);
b_A_x = x_vector(10);
b_A_y = x_vector(11);
b_A_z = x_vector(12);
b_p = x_vector(13);
b_q = x_vector(14);
b_r = x_vector(15);
V_wxE = x_vector(16);
V_wyE = x_vector(17);
V_wzE = x_vector(18);
A_x_m = c_m_vector(1);
A_y_m = c_m_vector(2);
A_z_m = c_m_vector(3);
p_m = c_m_vector(4);
q_m = c_m_vector(5);
r_m = c_m_vector(6);
H = [1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
end
