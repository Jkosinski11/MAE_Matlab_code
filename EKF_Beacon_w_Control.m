clear
clc
close all

function  Node_map = func1(letter)
Ts = 0.2;  % Time Step
N = 1000; % Total time steps

%% noise assumptions
sigma_V=0.2;
sigma_omega=0.1;
sigma_range=0.5;
sigma_phi=0.2;

%% Initial States
 X_Inital=[1 0 pi/2]'; % [X, Y, Heading]
% X_Inital=[5 -5 pi/2]'
%X_Inital=[10 0 pi]'

%% Initialize matrixes
Time=zeros(N,1);                                 % The time vector
V_T = zeros(N,1);                                % The true velocity input to the robot
omega_T = zeros(N,1);                            % The true anuglar rate input to the robot
V_I = zeros(N,1);                                % The noise-corrupted velocity input to the robot localization algorithm
omega_I = zeros(N,1);                            % The noise-corrupted anuglar rate input to the robot localization algorithm
x_T = zeros(N,1);                                % The truth for x position
y_T = zeros(N,1);                                % The truth for y position
theta_T = zeros(N,1);                            % The truth for robot heading
range_T = zeros(N,1);                            % The true distance between the beacon and the robot
phi_T = zeros(N,1);                              % The true bearing angle of the beacon with respect to the robot
range_M = zeros(N,1);                            % The measured distance between the beacon and the robot
phi_M = zeros(N,1);                              % The measured bearing angle of the beacon with respect to the robot
x_DR = zeros(N,1);                               % The Dead Reckoning estimated x-position
y_DR = zeros(N,1);                               % The Dead Reckoning y-position
theta_DR = zeros(N,1);                           % The Dead Reckoning robot heading
x_E = zeros(N,1);                                % The EKF estimated x-position
y_E = zeros(N,1);                                % The EKF estimated y-position
theta_E = zeros(N,1);                            % The EKF estimated robot heading
position_error_E=zeros(N,1);                     % distance from the truth (EKF)
position_error_DR=zeros(N,1);                    % distance from the truth (Dead reckoning)

goal_points = letter'; 
cur_goal_point = 0;
goal_threshold = 0.3;

% Gains for P controller
linear_P_gain = 1;
angular_P_gain = 5;

% Limits on velocities
linear_vel_lim = [0, 2];
angular_vel_lim = [-3, 3];


%% Initialization for the EKF
Q=Ts^2*[sigma_V^2   0               0
        0           sigma_V^2       0
        0           0               sigma_omega^2];                % Process noise covariance      
R= [sigma_range^2   0
    0               sigma_phi^2];                            % Measurement noise variance
I=eye(3);                                       % Identity matrix
X_P=X_Inital;                                % Initial estimate of the states (posterior)
X_P_DR=X_P;   
X_A=X_P;                                        % Initial estimate of the states (a prior)
x_DR(1)=X_A(1);
y_DR(1)=X_A(2);
theta_DR(1)=X_A(3);
x_E(1)=X_P(1);
y_E(1)=X_P(2);
theta_E(1)=X_P(3);
P_P=eye(3);                                     % Initial estimated error covariance

%% Main Simulation
% Each Time Step:
% Robot side:
% Get measurements from truth data and apply error
% Use DR to estimate state
% Use EKF to estimate state
% Calculate command linear and angular velocities using P controller
% If close enough to goal point, advance 1
% Simulation Side:
% Apply command velocities with kinematics

% Initialize Truth Data at first time step:
x_T(1) = X_Inital(1);
y_T(1) = X_Inital(2);
theta_T = X_Inital(3);

for k = 1:N
    %% Get measurement data from truth:
    range_T(k) = sqrt(x_T(k)^2 + y_T(k)^2);
    range_M(k) = range_T(k) + normrnd(0,sigma_range);    % create the sensor measurements
    phi_T(k)=atan2(y_T(k), x_T(k));
    phi_M(k)=phi_T(k) + normrnd(0,sigma_phi);
    
    %% Use EKF to estimate State:
    % Step #1: Predict states
    X_A(1)=X_P(1)+V_I(k)*cos(X_P(3))*Ts;        % Calculate the a priori (predicted) states 
    X_A(2)=X_P(2)+V_I(k)*sin(X_P(3))*Ts;        % Calculate the a priori (predicted) states 
    X_A(3)=X_P(3)+omega_I(k)*Ts;

    % Calculate the Jacobians for F and H
    F=[ 1       0           -Ts*V_I(k)*sin(X_A(3))
        0       1           Ts*V_I(k)*cos(X_A(3))
        0       0           1];
    R_Squre=X_A(1)^2+X_A(2)^2;                  % Temp variable that's equal to range squared
    H= [X_A(1)/sqrt(R_Squre)    X_A(2)/sqrt(R_Squre)    0
        -X_A(2)/R_Squre         X_A(1)/R_Squre          0];
    % Step #2: Predict the error covariance     
    P_A=F*P_P*F'+Q;                             % calcuate the a priori error covariance
    % Step #3: Calcuated the Kalman gain K
    K_k=P_A*H'*inv(H*P_A*H'+R);                 % Calcuate the Kalman Gain
    % Step #4: Update states:
    Residual=[range_M(k)-sqrt(R_Squre); wrapToPi(phi_M(k)-atan2(X_A(2),X_A(1)))]; % Calculate the residual
    X_P=X_A+K_k*Residual;  % calcuate the posterior states
    % Step #5: Update error covariance
    P_P=(I-K_k*H)*P_A;                          % calcuate the posterior error covariance
    % Get the EKF estimated states for plotting
    x_E(k)=X_P(1);
    y_E(k)=X_P(2);
    theta_E(k)=X_P(3);
    Theta_Sigma(k)=sqrt(P_P(3,3));
    position_error_E(k)=sqrt((x_E(k)-x_T(k))^2+(y_E(k)-y_T(k))^2);

    %% Use P Controller to command robot

    % Find distance error and heading error
    vect_to_goal = goal_points(:, cur_goal_point + 1) - [x_E(k); y_E(k)];
    dist_to_goal = norm(vect_to_goal);
    heading_to_goal = atan2(vect_to_goal(2), vect_to_goal(1));
    heading_error = heading_to_goal - mod(theta_E(k), 2*pi);
    if heading_error > pi
        heading_error = heading_error - 2*pi;
    end
    if heading_error < -pi
        heading_error = heading_error + 2*pi;
    end

    % Apply distance and heading errors with proportial gains:
    linear_vel_command = dist_to_goal * linear_P_gain;
    angular_vel_command = .98 *heading_error * angular_P_gain;

    % Determine if can advance to next waypoints
    if dist_to_goal < goal_threshold
        cur_goal_point = mod(cur_goal_point + 1, size(goal_points, 2));
    end

    %% Apply command velocities with kinematics:
    % Limit commanded velocities
    V_T(k+1) = min(linear_vel_lim(2), max(linear_vel_lim(1), linear_vel_command));
    omega_T(k+1) = min(angular_vel_lim(2), max(angular_vel_lim(1), angular_vel_command));
    V_I(k+1) = V_T(k+1) + normrnd(0,sigma_V);            % added noised to speed
    omega_I(k+1) = omega_T(k+1) + normrnd(0,sigma_omega);   % added noised to angular rate input

    % Propogate position
    x_T(k+1)=x_T(k)+V_T(k)*cos(theta_T(k))*Ts;
    y_T(k+1)=y_T(k)+V_T(k)*sin(theta_T(k))*Ts;
    theta_T(k+1)=theta_T(k)+omega_T(k)*Ts;

    %% Do some plotting
    plot(goal_points(1, :), goal_points(2, :), 'rs');
    hold on
    plot(0, 0, 'b*');
%    plot(x_T(k), y_T(k), 'bo', MarkerSize=20);
    plot(x_T(k), y_T(k), 'bo');
    plot(x_T(1:k), y_T(1:k), 'b');
    plot(x_E(1:k), y_E(1:k), 'r');
    hold off
    legend('Goal Points', 'Beacon', 'True Robot Position', 'True Robot Path', 'EKF Robot Path');
    grid minor
    xlabel('robot x position (m)');
    ylabel('robot y position (m)');
    title('Robot Navigation Using a Beacon')
    xlim([-18, 200]);
    ylim([-18, 200]);
    drawnow
    pause(0.01);
end
end

function data = readLetterCSV(letter)
   
    
    % Validate input
    if ~ischar(letter) || length(letter) ~= 1 || ~ismember(letter, 'A':'Z')
        error('Input must be a letter (A-Z).');
    end
    
    % Construct the filename
    filename = strcat('Documents\CPE412FinalProject\',letter, '.csv');
    
    % Check if the file exists
    if isfile(filename)
        % Read the CSV file
        data = readtable(filename);
        disp(['Successfully read file: ', filename]);
    else
        % Throw an error if the file does not exist
        error('File %s does not exist.', filename);
    end
end


prompt= 'input a String ';
string = upper(input(prompt, 's'));


for i = 1: length(string)
disp(string(i));
letter_Matrix = table2array(readLetterCSV(string(i)));
func1(letter_Matrix);
end