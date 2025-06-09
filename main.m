%%Extract the data
data = readmatrix('A1-007-DST-US06-FUDS-40-20120822.xlsx', 'Sheet', 'Sheet1');

%%Initialization
data = [data(:,1) data(:,4) data(:,5)]; %time = 1 Ut = 2, It = 3
time = data(:,1);
x = [data(:,3) data(:,2)]';
save('battery_data.mat', 'time', 'x');

%% Load Data
clear;
clear functions;
load('battery_data.mat');
x(2,:) = -x(2,:); % Invert current direction to match battery model (negative for discharge)

%% Coulomb Counting for True SOC and OCV
% Battery parameters
Qn = 1.1 * 3600; % 1100 mAh = 1.1 Ah = 3960 As (coulombs)
initial_SOC = ((x(1,1)) - 2.7)/(3.6 - 2.7); % Initial SOC (Assume related to nominal voltage)
coulombic_efficiency = 1;

% Initialize true SOC array
true_SOC = zeros(1, length(time));
true_SOC(1) = initial_SOC;

% Coulomb counting integration
for i = 2:length(time)
    dt = time(i) - time(i-1);
    current = x(2, i); % Current at time i
    
    % Coulomb counting: SOC(k) = SOC(k-1) - (dt * eta * I(k)) / Qn
    % Note: Positive current is discharge, negative is charge
    true_SOC(i) = true_SOC(i-1) - (dt * coulombic_efficiency * current) / Qn;
    
    % Bound SOC between 0 and 1
    true_SOC(i) = max(0, min(1, true_SOC(i)));
end

%% Plot Terminal Voltage and Current
figure('Name', 'Battery Input Data', 'Position', [100 100 1200 600]);
    
subplot(3,1,1);
plot(time, x(1,:), 'LineWidth', 1.5); 
title('Terminal Voltage', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
grid on;

subplot(3,1,2);
plot(time, x(2,:), 'LineWidth', 1.5);
title('Terminal Current', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Current [A]', 'FontSize', 12);
grid on;

subplot(3,1,3);
plot(time, true_SOC, 'LineWidth', 1.5);
title('Battery SOC', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('SOC', 'FontSize', 12);
grid on;

%% Initialize Parameters
error_MIUKF = zeros(1,22);
T = 5;
update_interval = 90;

N = length(time);
theta_history = zeros(6, N);
lambda_history = zeros(1, N);
electParams_history = zeros(5, N);
Ut_est_history = zeros(1,N);
SOC_est_history = zeros(1,N);
Uocv_est_history = zeros(1,N);

[theta_out, lambda] = VFFRLS(x(:,1), error_MIUKF);
[x_pred, error_MIUKF, electParams] = MIUKF(x(:,1), theta_out, T);

%logging
theta_history(:,1) = theta_out;
lambda_history(1) = lambda;
electParams_history(:,1) = electParams;
Ut_est_history(1) = x_pred(1);
SOC_est_history(1) = x_pred(3);
Uocv_est_history(1) = theta_out(1)/(1 + theta_out(2) + theta_out(3));

for i = 2:N
    T = time(i) - time(i-1);
    if mod(i, update_interval) == 0
        [theta_out, lambda] = VFFRLS(x_pred, error_MIUKF);
    end
    [x_pred, error_MIUKF, electParams] = MIUKF(x(:,i), theta_out, T);
    
    %logging
    theta_history(:,i) = theta_out;
    lambda_history(i) = lambda;
    electParams_history(:,i) = electParams;
    Ut_est_history(i) = x_pred(1);
    SOC_est_history(i) = x_pred(3);
    Uocv_est_history(i) = theta_out(1)/(1 + theta_out(2) + theta_out(3));
end

%% Plotting Ut, Uocv, and SOC Estimates
figure('Name', 'Battery State Estimates', 'Position', [100 50 1200 700]);
    
subplot(3,1,1);
plot(time, Ut_est_history, 'b-', 'LineWidth', 1.5);
title('Terminal Voltage Estimate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
grid on;

subplot(3,1,2);
plot(time, Uocv_est_history, 'r-', 'LineWidth', 1.5);
title('Open Circuit Voltage Estimate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
grid on;

subplot(3,1,3);
plot(time, SOC_est_history, 'g-', 'LineWidth', 1.5);
title('State of Charge (SOC) Estimate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('SOC', 'FontSize', 12);
grid on;

%% Plotting Lambda History
figure('Name', 'Variable Forgetting Factor (VFF) History', 'Position', [100 50 1200 700]);
plot(time, lambda_history, 'k-', 'LineWidth', 1.5);
title('Variable Forgetting Factor (VFF) History', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Lambda', 'FontSize', 12);
grid on;

%% Plotting Electrical Parameters
figure('Name', 'Electrical Parameters', 'Position', [100 50 1200 700]);
    
param_names = {'R0', 'R1', 'C1', 'R2', 'C2'};
units = {'Ω', 'Ω', 'F', 'Ω', 'F'};
colors = {'b', 'r', 'g', 'm', 'c'};

for i = 1:5
    subplot(3,2,i);
    plot(time, electParams_history(i,:), colors{i}, 'LineWidth', 1.5);
    title(param_names{i}, 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 12);
    ylabel([param_names{i} ' [' units{i} ']'], 'FontSize', 12);
    grid on;
end

%% Plotting VFFRLS Parameters
figure('Name', 'VFFRLS Parameters', 'Position', [100 50 1200 700]);

param_names_vffrls = {'(1 + a1 + a2) * Uocv','a1', 'a2', 'b0', 'b1', 'b2'};
units_vffrls = {'', '', '', '', '', ''};
colors_vffrls = {'b', 'r', 'g', 'm', 'c', 'y'};
for i = 1:6
    subplot(3,2,i);
    plot(time, theta_history(i,:), colors_vffrls{i}, 'LineWidth', 1.5);
    title(param_names_vffrls{i}, 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 12);
    ylabel([param_names_vffrls{i} ' [' units_vffrls{i} ']'], 'FontSize', 12);
    grid on;
end

%% Plotting Clean SOC vs OCV with SOC=1 on Left

[sortedSOC, sortIndex] = sort(SOC_est_history);
sortedOCV = Uocv_est_history(sortIndex);

smoothedOCV = movmean(sortedOCV,1400);

% Create new figure with SOC=1 on the left
figure('Name', 'SOC vs Open Circuit Voltage (LiFePO4)', 'Position', [100 50 1200 700]);

% Plot with SOC=1 on the left by using 1-SOC
plot(1-sortedSOC, smoothedOCV, 'b-', 'LineWidth', 2.5);

% Improve graph appearance
title('SOC vs Open Circuit Voltage (LiFePO4)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('State of Charge (SOC)', 'FontSize', 14);
ylabel('Open Circuit Voltage (OCV) [V]', 'FontSize', 14);
grid on;
box on;

% Create better x-axis ticks for clearer reading
set(gca, 'FontSize', 12, 'XTick', 0:0.1:1);
xticklabels({'1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'});
ylim([3.1 3.6]);

%% Comparing Ut and SOC Estimates
figure('Name', 'Ut and SOC Estimates vs Real', 'Position', [100 50 1200 700]);
    
subplot(3,1,1);
plot(time, x(1,:), 'r-', 'LineWidth', 1.5);
hold on;
plot(time, Ut_est_history, 'b-', 'LineWidth', 1.5);

title('Terminal Voltage Estimate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
legend('Measured V_t', 'Estimated V_t', 'Location', 'best', 'FontSize', 10);
grid on;

subplot(3,1,2);
plot(time, SOC_est_history, 'g-', 'LineWidth', 1.5);
hold on;
plot(time, true_SOC, 'y-', 'LineWidth', 1.5);
title('State of Charge (SOC) Estimate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('SOC', 'FontSize', 12);
legend('Estimated SOC', 'True SOC', 'Location', 'best', 'FontSize', 10);
grid on;

subplot(3,1,3);
plot(time, abs(true_SOC - SOC_est_history), 'g-', 'LineWidth', 1.5);
title('State of Charge (SOC) Estimate Error', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Error', 'FontSize', 12);
grid on;

%% Save Graphs
% Create folder for saving if it doesn't exist
savePath = 'graphs';
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

% Get all valid figure handles
figHandles = findobj('Type', 'figure');

% Loop through each figure and save
for i = 1:length(figHandles)
    fig = figHandles(i);
    
    % Get figure name or number
    if isempty(fig.Name)
        figureName = sprintf('Figure%d', fig.Number);
    else
        figureName = strrep(fig.Name, ' ', '_');
        figureName = strrep(figureName, ':', '');
        figureName = strrep(figureName, '/', '_');
    end
    
    % Save as PNG (more robust method)
    pngFile = fullfile(savePath, [figureName '.png']);
    
    % Check if figure is valid
    if isvalid(fig)
        % Use print instead of saveas for better reliability
        print(fig, pngFile, '-dpng', '-r300');
        fprintf('Saved: %s\n', pngFile);
    else
        fprintf('Skipping invalid figure: %s\n', figureName);
    end
end

fprintf('All figures saved to "%s" folder\n', savePath);