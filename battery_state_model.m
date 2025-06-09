function pred = battery_state_model (u, x, T, electParams)
    %Extract battery parameters
    R1 = electParams(2);
    C1 = electParams(3);
    R2 = electParams(4);
    C2 = electParams(5);

    % Extract States
    U1 = x(1); % Voltage at the first capacitor
    U2 = x(2); % Voltage at the second capacitor
    SOC = x(3); % State of Charge
    I = u; % Current input

    %Constant parameters
    Qn = 1.1 * 3600; % Nominal capacity in As (example value, adjust as needed)
    n = 1; %coulomb efficiency (example value, adjust as needed)

    % Update the battery state based on the input current and parameters
    tau1 = R1 * C1;
    tau2 = R2 * C2;

    % Forward Euler Modelling
    if T < tau1
        predU1 = (1 - T/tau1) * U1 + (T / C1) * I;
    else
        predU1 = R1 * I; % If T exceeds tau1, assume steady state
    end
    
    if T < tau2
        predU2 = (1 - T/tau2) * U2 + (T / C2) * I;
    else
        predU2 = R2 * I; % If T exceeds tau2, assume steady state
    end

    predSOC = SOC - (T * n / Qn) * I; % Update SOC based on current
    predSOC = max(0, min(1, predSOC)); % Ensure SOC is within [0, 1]

    pred = [predU1 predU2 predSOC]; % Predicted state vector
end