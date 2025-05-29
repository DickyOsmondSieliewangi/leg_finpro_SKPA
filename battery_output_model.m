function y = battery_output_model(x, u, electParams, theta_out)
    % Extract states
    U1 = x(1); % Voltage at the first capacitor
    U2 = x(2); % Voltage at the second capacitor
    I = u; % Current input

    % Extract electrical parameters
    R0 = electParams(1);
    a1 = theta_out(2);
    a2 = theta_out(3);
    Uocv = theta_out(1) / (1 + a1 + a2); % Open circuit voltage

    % Calculate the output voltage based on the battery model
    y = Uocv - U1 - U2 - I * R0; % Example output calculation
end