function electParams = Levenberg_Marquardt(theta_out, T)
    %Extract parameters from theta_out
    a1 = theta_out(2);
    a2 = theta_out(3);
    b0 = theta_out(4);
    b1 = theta_out(5);
    b2 = theta_out(6);

    %Calculate R0, R1, R2, C1, C2
    
    %a1 = (T - 2* C2 * R2) / (T + 2 * C2 * R2) + (T - 2* C1 * R1) / (T + 2 * C1 * R1);
    %a2 = ((T - 2* C1 * R1) * (T - 2* C2 * R2)) / ((T + 2 * C1 * R1) * (T + 2 * C2 * R2));
    %bo = T * (R1/(T + 2 * C1 * R1) + R2/(T + 2 * C2 * R2)) + R0;
    %b1 = (2 * T ^ 2 * (R1 + R2) + R0 * (2 * T ^ 2 - 8 * C1 * C2 * R1 * R2)) / ((T + 2 * C1 * R1) * (T + 2 * C2 * R2));
    %b2 = (T * (R1 * (T - 2 * C2 * R2) + R2 * (T - 2 * C1 * R1)) + R0 * (T - 2 * C2 * R2) * (T - 2 * C1 * R1)) / ((T + 2 * C1 * R1) * (T + 2 * C2 * R2));

    % Initial guess for electrical parameters (from paper)
    persistent prev_electParams;
    if isempty(prev_electParams)
        prev_electParams = [0.002, 0.0012, 7.32e4, 0.0011, 4.49e4];
    end
    x0 = prev_electParams; % [R0, R1, C1, R2, C2]
    
    % Define the nonlinear system of equations
    equations = @(x) residual_function(x, a1, a2, b0, b1, b2, T);
    
    % Set optimization options
    options = optimoptions('lsqnonlin', ...
        'Algorithm', 'levenberg-marquardt', ...
        'Display', 'off', ...
        'TolFun', 1e-5, ...
        'TolX', 1e-5, ...
        'MaxIter', 1000);
    
    % Set bounds (all parameters should be positive)
    lb = [0, 0, 0, 0, 0];  % Lower bounds
    ub = [];  % Upper bounds
    
    % Solve the nonlinear system
    [x_opt, ~, ~, exitflag] = lsqnonlin(equations, x0, lb, ub, options);
    
    % MATLAB automatically:
    % 1. Calculates Jacobian J
    % 2. Solves (J'J + μI)δ = -J'r
    % 3. Updates parameters: [R0, R1, C1, R2, C2] = x + δ
    % 4. Adjusts μ based on improvement

    if exitflag <= 0
        warning('Levenberg-Marquardt failed to converge');
        electParams = prev_electParams;
        return % Use previous parameters if failed
    end
    
    % Extract results
    R0 = x_opt(1);
    R1 = x_opt(2);
    C1 = x_opt(3);
    R2 = x_opt(4);
    C2 = x_opt(5);

    electParams = [R0 R1 C1 R2 C2];
end

function residuals = residual_function(x, a1_target, a2_target, b0_target, b1_target, b2_target, T)
    % Extract current parameter estimates
    R0 = x(1);
    R1 = x(2);
    C1 = x(3);
    R2 = x(4);
    C2 = x(5);
    
    % Calculate what a1, a2, b0, b1, b2 should be with current parameters
    % Using your forward equations:
    tau1 = C1 * R1; % Time constant for first RC pair
    tau2 = C2 * R2; % Time constant for second RC pair
    
    a1_calc = (T - 2*tau2) / (T + 2*tau2) + (T - 2*tau1) / (T + 2*tau1);
    
    a2_calc = ((T - 2*tau1) * (T - 2*tau2)) / ((T + 2*tau1) * (T + 2*tau2));
    
    b0_calc = T * (R1/(T + 2*tau1) + R2/(T + 2*tau2)) + R0;
    
    b1_calc = (2*T^2*(R1 + R2) + R0*(2*T^2 - 8*tau1*tau2)) / ...
              ((T + 2*tau1) * (T + 2*tau2));
    
    b2_calc = (T*(R1*(T - 2*tau2) + R2*(T - 2*tau1)) + ...
               R0*(T - 2*tau2)*(T - 2*tau1)) / ...
              ((T + 2*tau1) * (T + 2*tau2));
    
    % Calculate residuals (difference between target and calculated values)
    residuals = [
        a1_calc - a1_target;
        a2_calc - a2_target;
        b0_calc - b0_target;
        b1_calc - b1_target;
        b2_calc - b2_target
    ];
end