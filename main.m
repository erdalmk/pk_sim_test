%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min
true_params.ke_Central = 3e-2;  % 1/min

%% Get Simulated Data
[time, u, y_conc, x_conc] = get_data(true_params);

%% Get Fitting Structure
fitObj = FittingStruct(time, u, y_conc, 'Exact');

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

%% Display Results
% Plot simulated response
figure('Units', 'normalized', 'OuterPosition', [0.1, 0.1, 0.4, 0.8])
subplot(3,1,1:2)
plot(time, y_conc, '.', 'MarkerSize',8)
hold on
plot(fitObj.time, sol_nlp.x, 'Color', .7*[0,1,0], 'LineWidth',2)
xlim(time([1, end]))

legend('data','fit')

ylabel('Concentration (\mu M)')
xlabel('Time (min)')

set(gca,'TickDir', 'out', 'box', 'off')

subplot(3,1,3)
stairs(time, u, 'm')
ylabel('Injection Rate (\mu mol/min)')

xlim(time([1, end]))
xlabel('Time (min)')

set(gca,'TickDir', 'out', 'box', 'off')

% Show estimated parameters against true parameters
T = table;
T.parameter = fields(true_params);
T.true = zeros(4, 1);
T.estimated = zeros(4, 1);
for i = 1:4
    T.true(i) = true_params.(T.parameter{i});
    T.estimated(i) = sol_nlp.(T.parameter{i});
end
disp(T)

