function [theta_out, lambda, y_pred] = VFFRLS(x, error_MIUKF)
    %From x calculate phi
    %x = [U1(k) U2(k) SOC(k)]
    %phi = [1 -Ut(k-1) -Ut(k-2) -I(k) -I(k-1) -I(k-2)]

    newUt = x(1);
    newI = x(2);

    % Initialize Ut and I if they are empty
    persistent Ut I;
    if isempty(Ut)
        Ut = zeros(1,3);
    end

    if isempty(I)
        I = zeros(1,3);
    end

    % Shift the values to make room for the new ones
    Ut = [newUt Ut(1:end-1)];
    I = [newI I(1:end-1)];

    phi = [1 -Ut(2) -Ut(3) -I(1) -I(2) -I(3)]';
    n = length(phi);

    %From error_MIUKF calculate lambda
    S = 22;
    rho = 200;
    lambda_min = 0.9;
    lambda_max = 0.995;

    sum_squared_error = error_MIUKF * error_MIUKF'; % Calculate the sum of squared errors

    % Calculate L(k+1) = -ρ * Σ(ei * ei^T) / S
    L = -rho * sum_squared_error / S;
    
    % Calculate variable forgetting factor
    lambda = lambda_min + (lambda_max - lambda_min) * (2^L);
    
    % Ensure lambda is within bounds
    lambda = max(lambda_min, min(lambda_max, lambda));

    % Initialize theta and P if they are empty
    persistent theta P;
    if isempty(theta)
        theta = zeros(n, 1);
    end

    if isempty(P)
        P = (1e2) * eye(n);
    end

    denominator = lambda + phi' * P * phi;

    % Calculate the Algorithm gain
    K = P * phi / denominator;

    %Calculate covariance matrix
    P = (P - K * phi' * P) / lambda;

    %Predict y
    y_pred = phi' * theta;

    %Check error
    error = newUt - y_pred;

    %Update theta
    theta_new = theta + K * error;

    theta = theta_new;

    %Output theta
    theta_out = theta;
end