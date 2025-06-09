function [x_pred, error_MIUKF, electParams] = MIUKF(x, theta_out, T)
    % MIUKF: Multi Innovation Unscented Kalman Filter

    % extract inputs
    y = x(1);
    u = x(2);

    %Use levenberg-marquardt to calculate electrical parameters
    electParams = Levenberg_Marquardt(theta_out, T);

    %Calculate Sampling Points
    persistent posterior_x posterior_P;
    if isempty(posterior_x)
        posterior_x = zeros(3, 1); % Initial state vector [U1, U2, SOC]
        posterior_x(3) = (2.8 - 2.7)/(3.6 - 2.7); % Initial SOC, can be adjusted based on prior knowledge
    end
    if isempty(posterior_P)
        posterior_P = eye(3) * 1e-3; % Initial covariance matrix
    end

    L = 3; % Number of state variables and sigma points
    alpha = 1e-3; % Scaling parameter can be 1e-3
    beta = 2.0; % UKF Params
    ki = 0; % UKF Params

    eta = alpha^2 * (L + ki) - L; % Scaling ratio

    weights_x = zeros(2 * L + 1, 1); % State weights
    weights_x(1) = eta / (eta + L); % First weight is for the central sigma point
    weights_x(2:end) = 1 / (2 * (eta + L)); % Weights for the other sigma points

    weights_P = zeros(2 * L + 1, 1); % Covariance weights
    weights_P(1) = eta / (eta + L) + (1 - alpha^2 + beta); % First weight is for the central sigma point
    weights_P(2:end) = 1 / (2 * (eta + L)); % Weights for the other sigma points

    % Generate sigma points
    sigmaPoints = zeros(2 * L + 1, L); % 3 state variables: U1 U2 and SOC
    sigmaPoints(1, :) = posterior_x'; % First sigma point is the mean
    
    % In MIUKF, before Cholesky decomposition:
    try
        sqrtFactor = chol((eta + L) * posterior_P, 'lower');
    catch
        warning('posterior_P not positive definite — using filter');
        posterior_P_mod = (posterior_P + posterior_P') / 2;                     % Symmetrize
        posterior_P_mod = posterior_P_mod + 1e-8 * eye(size(posterior_P_mod));       % Regularize
        sqrtFactor = chol((eta + L) * posterior_P_mod, 'lower');
    end
    
    % Add each column of sqrtFactor to get positive perturbations
    for i = 1:L
        sigmaPoints(1 + i, :) = posterior_x' + sqrtFactor(:, i)'; % Add i-th column
    end
    % Subtract each column of sqrtFactor to get negative perturbations  
    for i = 1:L
        sigmaPoints(1 + L + i, :) = posterior_x' - sqrtFactor(:, i)'; % Subtract i-th column
    end

    % Predict sigma points
    % For each sigma point, we will calculate the predicted state
    predictedSigmaPoints = zeros(size(sigmaPoints)); % Initialize predicted sigma points
    for i = 1:size(sigmaPoints, 1)
        predictedSigmaPoints(i, :) = battery_state_model(u, sigmaPoints(i, :), T, electParams);
    end

    %Update apriori state
    % Calculate the a priori for state and covariance
    a_priori_x = predictedSigmaPoints' * weights_x; % Weighted mean of predicted sigma points
    a_priori_P = zeros(size(posterior_P)); % Initialize a priori covariance
    deviation = predictedSigmaPoints - a_priori_x'; % Deviation from the mean
    for i = 1:2*L+1
        diff = deviation(i, :)'; % [3×1] deviation vector
        a_priori_P = a_priori_P + weights_P(i) * (diff * diff'); % [3×3] += scalar * [3×3]
    end

    Qk = 1e-8 * eye(3); % Process noise covariance matrix, assuming small process noise
    a_priori_P = a_priori_P + Qk; % Add the previous covariance

    % Calculate the observation a priori value, 
    % the prediction value of observation variance, 
    % and estimate the covariance difference and gain

    % Assuming y is the observation vector, we will calculate the predicted observation
    y_sigma = zeros(size(sigmaPoints, 1), 1); % Initialize Outputs
    for i = 1:size(sigmaPoints, 1)
        %Measurement model for each sigma point
        %y = Uocv - U1 - U2 - I * R0
        y_sigma(i) = battery_output_model(predictedSigmaPoints(i, :)', u, electParams, theta_out);
    end

    y_pred = y_sigma' * weights_x; % Weighted mean of predicted outputs
    
    obs_P = 0; % Initialize observation covariance Pyy

    for i = 1:2*L+1
        diff_y = y_sigma(i) - y_pred; % Measurement deviation
        obs_P = obs_P + weights_P(i) * diff_y^2; % Accumulate weighted variance
    end

    % Add measurement noise
    R = 1e-8; % Measurement noise variance (terminal voltage noise)
    obs_P = obs_P + R;

    diff_P = zeros(L, 1); % Initialize observation covariance Pxy % [3×1] for 3 states, 1 measurement
    for i = 1:2*L+1
        diff_x = predictedSigmaPoints(i, :)' - a_priori_x; % [3×1] state deviation
        diff_y = y_sigma(i) - y_pred; % Scalar measurement deviation
        diff_P = diff_P + weights_P(i) * diff_x * diff_y; % [3×1] += scalar * [3×1] * scalar
    end

    % Calculate the Kalman gain
    K = diff_P / obs_P; % [3×1] gain vector

    %Calculate the multi information error  
    % Assuming y is the voltage at the terminal
    M = 22;
    a = 0.5;

    persistent error_history gain_history;
    if isempty(error_history)
        error_history = zeros(1, M); % Initialize error history
    end
    if isempty(gain_history)
        gain_history = zeros(length(K), M); % Initialize gain history
    end

    currentError = y - y_pred; % Innovation or measurement residual
    error_history = [currentError error_history(1:end-1)]; % Shift the error history / innovation
    gain_history = [K gain_history(:,1:end-1)]; % Shift the gain history

    gamma = zeros(1, M); % Initialize gamma vector
    gamma(1) = 1;
    for i = 2:M
        gamma(i) = a / (M - 1);
    end
    %matrix by  measured value, and M is the information length

    errorMI = 0;
    for i = 1:M
        errorMI = errorMI + gamma(i) * gain_history(:, i) * error_history(i);
        %[3×1]    = [3×1] + scalar   * [3×1]           * scalar
    end

    %Update posterior state
    posterior_x = a_priori_x + errorMI;
    posterior_x(3) = max(0, min(1, posterior_x(3)));
    posterior_P = a_priori_P - K * obs_P * K';

    %Outputs
    Ut = battery_output_model(posterior_x, u, electParams, theta_out);
    I = u;
    SOC = posterior_x(3);
    x_pred = [Ut I SOC];
    error_MIUKF = error_history;
end