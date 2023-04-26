%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min
true_params.ke_Central = 3e-2;  % 1/min

%% Get Simulated Data
[time, u, y_conc, x_conc] = get_data(true_params);

%% Get Fitting Structure
fitObj = FittingStruct(time, u, y_conc);

%% Run Solver
% Initial Values
V1_init = 1e-2;
V2_init = 1e-3;
CL12_init = 1e-5;
kE_init = 1e-4;

% initial rates
init_vals = struct;
init_vals.kU = 1/V1_init;
init_vals.kE = kE_init;
init_vals.k12 = CL12_init/V2_init;
init_vals.k21 = CL12_init/V1_init;

init_vals.invsig1 = 1e1;
init_vals.invsig2 = 1e1;

sol_nlp = fitObj.fit_nlp(init_vals);

%% Plot simulated response
figure
plot(time, y_conc, '.')
hold on
set(gca,'ColorOrderIndex',1)
plot(time, x_conc)
plot(fitObj.time, sol_nlp.x)
ylabel('Concentration (\mu M)')

yyaxis right
stairs(time, u, 'r')
ylabel('Injection Rate (\mu mol/min)')

xlim(time([1, end]))
xlabel('Time (min)')

ax = gca;
ax.YAxis(2).Color = 'r';
set(ax,'TickDir', 'out', 'box', 'off')