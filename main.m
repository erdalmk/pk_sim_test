%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min
true_params.ke_Central = 3e-2;  % 1/min

%% Get Simulated Data
[time, y_conc, x_conc] = get_data(true_params);

%% Plot simulated response
figure
plot(time, y_conc, '.')
hold on
set(gca,'ColorOrderIndex',1)
plot(time, x_conc)
xlim(time([1, end]))
ylim([-20, 200])
set(gca,'TickDir', 'out', 'box', 'off')