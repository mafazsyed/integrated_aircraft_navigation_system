% Template for Kalman Filter (KF), Extended Kalman Filter (EKF) and Iterated Extended Kalman Filter (IEKF)

% clear; clear all; close all;

% Read data files
load('dataTask1.mat')
z_k=d_k;
u_k = c_k;
% Initialization
Ts=dt;     %time step (already provided by the data)
N=length(t);      %total number of steps
% standard deviation of input measurement noise w
stdw = [0.01 0.01 0.01 0.01*pi/180 0.01*pi/180 0.01*pi/180];
% standard deviation of output measurement noise v
stdv = [5 5 10 0.1 0.1 0.1 0.1*pi/180 0.1*pi/180 0.1*pi/180 0.1 0.1*pi/180 0.1*pi/180];
Ex_0=[z_k(1,1) z_k(1,2) z_k(1,3) 95 0 0 z_k(1,7) z_k(1,8) z_k(1,9) 0 0 0];      %expected value of x_0
stdx_0 = [1 1 1 50 50 50 1 1 1 50 50 50];  %standard deviation of x_0


xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
P_km1_km1 = diag(stdx_0.^2);  % P(0|0) = P(0)

Q = diag(stdw.^2);
R = diag(stdv.^2);

n = length(xhat_km1_km1); % n: state dimension
m = size(u_k, 2);     % m: observation dimension
p = size(z_k, 2);     % m: observation dimension
u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices

% Preallocate storage
stdx_cor  = zeros(N, n);  % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
x_cor     = zeros(N, n);  % \hat{x}(k-1|k-1), previous estimation
K       = cell(N, 1);   % K(k) Kalman Gain
innov     = zeros(N, p);  % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);

for k=1:N
    % Step 1: Prediction
    [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
    xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)

    % Step 2: Covariance matrix of state prediction error / Minimum prediction MSE
    [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
    P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)

    % Step 3: Kalman Gain
    H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
    Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
    K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)
    
    % Step 4: Measurement Update (Correction)
    z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction)
    xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k'; % x(k|k) (correction)

    % Step 5: Correction for Covariance matrix of state Estimate error /
    % Minimum MSE
    I_KH = eye(n) - K_k * H_k;
    P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)

    % Save data: State estimate and std dev
    stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
    x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
    K{k,1} = K_k; % K(k) (gain)
    innov(k,:)= z_k(k,:) - z_k_km1;
    
    % Recursive step
    xhat_km1_km1 = xhat_k_k; 
    P_km1_km1 = P_k_k;

end


%% Plots
x_cor_task1 = x_cor;

figure
subplot(6,2,1)
plot(t(1:N),x_cor(1:N,1),'b','LineWidth',2)
title('x_{E}')
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(6,2,2)
plot(t(1:N),x_cor(1:N,2),'b','LineWidth',2)
title('y_{E}')
xlabel('Time (s)')
ylabel('Position (m)')
grid on
% legend('raw measurements','original','estimation')

subplot(6,2,3)
plot(t(1:N),x_cor(1:N,3),'b','LineWidth',2)
title('z_{E}')
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(6,2,4)
plot(t(1:N),x_cor(1:N,4),'b','LineWidth',2)
title('u')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

subplot(6,2,5)
plot(t(1:N),x_cor(1:N,5),'b','LineWidth',2)
title('v')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

subplot(6,2,6)
plot(t(1:N),x_cor(1:N,6),'b','LineWidth',2)
title('w')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

subplot(6,2,7)
plot(t(1:N),x_cor(1:N,7),'b','LineWidth',2)
title('{\phi}')
xlabel('Time (s)')
ylabel('Angle (rad)')
grid on

subplot(6,2,8)
plot(t(1:N),x_cor(1:N,8),'b','LineWidth',2)
title('{\theta}')
xlabel('Time (s)')
ylabel('Angle (rad)')
grid on

subplot(6,2,9)
plot(t(1:N),x_cor(1:N,9),'b','LineWidth',2)
title('{\psi}')
xlabel('Time (s)')
ylabel('Angle (rad)')
grid on

subplot(6,2,10)
plot(t(1:N),x_cor(1:N,10),'b','LineWidth',2)
title('V_{wxE}')
xlabel('Time (s)')
ylabel('Wind Velocity (m/s)')
grid on

subplot(6,2,11)
plot(t(1:N),x_cor(1:N,11),'b','LineWidth',2)
title('V_{wyE}')
xlabel('Time (s)')
ylabel('Wind Velocity (m/s)')
grid on

subplot(6,2,12)
plot(t(1:N),x_cor(1:N,12),'b','LineWidth',2)
title('V_{wzE}')
xlabel('Time (s)')
ylabel('Wind Velocity (m/s)')
grid on

%% Functions

function x_state_dot = funcf(xi,u,t)
    
    % Extract values
    A_x_m = u(1);
    A_y_m = u(2);
    A_z_m = u(3);
    p_m = u(4);
    q_m = u(5);
    r_m= u(6);

    x = xi(1) ;
    y = xi(2) ;
    z = xi(3);
    u = xi(4);
    v = xi(5);
    w = xi(6);
    phi = xi(7);
    theta = xi(8);
    psi = xi(9);
    v_wxE = xi(10);
    v_wyE = xi(11);
    v_wzE = xi(12);
    
    
    
    
    %Bias Free Kinematic Model (excluding noises)
    x_dot =v_wxE - sin(psi)*(v*cos(phi) - w*sin(phi)) + cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta));
    y_dot =v_wyE + sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) + cos(psi)*(v*cos(phi) - w*sin(phi));
    z_dot =v_wzE - u*sin(theta) + cos(theta)*(w*cos(phi) + v*sin(phi));
    u_dot =A_x_m - (981*sin(theta))/100 - q_m*w + r_m*v;
    v_dot =A_y_m + (981*cos(theta)*sin(phi))/100 + p_m*w - r_m*u;
    w_dot =A_z_m + (981*cos(phi)*cos(theta))/100 - p_m*v + q_m*u;
    phi_dot =p_m + r_m*cos(phi)*tan(theta) + q_m*sin(phi)*tan(theta);
    theta_dot =q_m*cos(phi) - r_m*sin(phi);
    psi_dot =(r_m*cos(phi))/cos(theta) + (q_m*sin(phi))/cos(theta);
    v_wxE_dot =0;
    v_wyE_dot =0;
    v_wzE_dot =0;

    x_state_dot = [x_dot; y_dot; z_dot; u_dot; v_dot; w_dot; phi_dot; theta_dot; psi_dot; v_wxE_dot; v_wyE_dot; v_wzE_dot];
end

function y = funch(xi,u,t)

    % Extract values
    A_x = u(1) ;
    A_y = u(2) ;
    A_z = u(3) ;
    p = u(4) ;
    q = u(5) ;
    r= u(6) ;
    x = xi(1) ;
    y = xi(2) ;
    z = xi(3);
    u = xi(4);
    v = xi(5);
    w = xi(6);
    phi = xi(7);
    theta = xi(8);
    psi = xi(9);
    v_wxE = xi(10);
    v_wyE = xi(11);
    v_wzE = xi(12);
    
    
    x_GPS=x;
    y_GPS=y;
    z_GPS=z;
    u_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    v_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    w_GPS=-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    phi_GPS=phi;
    theta_GPS=theta;
    psi_GPS=psi;
    V=sqrt(u^2+v^2+w^2);
    alpha=atan(w/u);
    beta=atan(v/sqrt(u^2+w^2));
    y = [x_GPS, y_GPS, z_GPS, u_GPS, v_GPS, w_GPS, phi_GPS, theta_GPS, psi_GPS, V, alpha, beta];

end

function [Phi,Gamma] = funcLinDisDyn(x_linpt,u_linpt,Ts)

    % Extract values
    x = x_linpt(1) ;
    y = x_linpt(2) ;
    z = x_linpt(3);
    u = x_linpt(4);
    v = x_linpt(5);
    w = x_linpt(6);
    phi = x_linpt(7);
    theta = x_linpt(8);
    psi = x_linpt(9);
    v_wxE = x_linpt(10);
    v_wyE = x_linpt(11);
    v_wzE = x_linpt(12);
    A_x_m = u_linpt(1) ;
    A_y_m = u_linpt(2) ;
    A_z_m = u_linpt(3) ;
    p_m = u_linpt(4) ;
    q_m = u_linpt(5) ;
    r_m= u_linpt(6) ;
   
   
        
    % Numerical evaluation of continuous - time dynamics
    A=[[0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                 -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0]
[0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                 -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0]
[0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1]
[0, 0, 0,                   0,                                              r_m,                                             -q_m,                                                                                  0,                                                           -(981*cos(theta))/100,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                -r_m,                                                0,                                              p_m,                                                      (981*cos(phi)*cos(theta))/100,                                                  -(981*sin(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                 q_m,                                             -p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                  -(981*cos(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                                  q_m*cos(phi)*tan(theta) - r_m*sin(phi)*tan(theta),               r_m*cos(phi)*(tan(theta)^2 + 1) + q_m*sin(phi)*(tan(theta)^2 + 1),                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                                                      - r_m*cos(phi) - q_m*sin(phi),                                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                              (q_m*cos(phi))/cos(theta) - (r_m*sin(phi))/cos(theta), (r_m*cos(phi)*sin(theta))/cos(theta)^2 + (q_m*sin(phi)*sin(theta))/cos(theta)^2,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0]];
 

    G = [[ 0,  0,  0,  0,                    0,                    0]
[ 0,  0,  0,  0,                    0,                    0]
[ 0,  0,  0,  0,                    0,                    0]
[-1,  0,  0,  0,                    w,                   -v]
[ 0, -1,  0, -w,                    0,                    u]
[ 0,  0, -1,  v,                   -u,                    0]
[ 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta)]
[ 0,  0,  0,  0,            -cos(phi),             sin(phi)]
[ 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta)]
[ 0,  0,  0,  0,                    0,                    0]
[ 0,  0,  0,  0,                    0,                    0]
[ 0,  0,  0,  0,                    0,                    0]];

    
    % Discretisation of dynamics
    [Phi , Gamma ]= c2d (A ,G , Ts ) ;

end

function H = funcLinDisObs(x_linpt,u_linpt,t)

    % Extract values
    x = x_linpt(1) ;
    y = x_linpt(2) ;
    z = x_linpt(3);
    u = x_linpt(4);
    v = x_linpt(5);
    w = x_linpt(6);
    phi = x_linpt(7);
    theta = x_linpt(8);
    psi = x_linpt(9);
    v_wxE = x_linpt(10);
    v_wyE = x_linpt(11);
    v_wzE = x_linpt(12);
    A_x = u_linpt(1) ;
    A_y = u_linpt(2) ;
    A_z = u_linpt(3) ;
    p = u_linpt(4) ;
    q = u_linpt(5) ;
    r= u_linpt(6) ;

    % Numerical evaluation
    H = [[1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0]
[0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0]
[0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1]
[0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0]
[0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]
[0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0]];

end
