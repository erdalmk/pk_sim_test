%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min
true_params.ke_Central = 3e-2;  % 1/min

%% Get Simulated Data
[time, u, y_conc, x_conc] = get_data(true_params);

%% Plot simulated response
figure
plot(time, y_conc, '.')
hold on
set(gca,'ColorOrderIndex',1)
plot(time, x_conc)
ylabel('Concentration (\mu M)')

yyaxis right
stairs(time, u, 'r')
ylabel('Injection Rate (\mu mol/min)')

xlim(time([1, end]))
xlabel('Time (min)')

ax = gca;
ax.YAxis(2).Color = 'r';
set(ax,'TickDir', 'out', 'box', 'off')