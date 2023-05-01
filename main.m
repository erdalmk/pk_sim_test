%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min

%% Set Simulation and Fitting Parameters
sim_config = struct;
sim_config.model_elim = 'enzymatic';
sim_config.uniform_sampling = true;
sim_config.sampling_time = 1; % ignored if uniform_sampling is false
sim_config.repeat_count = 3;
sim_config.SNR = 30;

disc_type = 'Exact';    % can be ForwardEuler or Exact
elim_type = 'varying';  % can be varying or constant

%% Set true parameters based on elimination type
if strcmp(sim_config.model_elim, 'enzymatic')
    true_params.Km_Central = 5e1;  % umole/L
    true_params.Vm_Central = 5e-2;  % umole/min
elseif strcmp(sim_config.model_elim, 'linear')
    true_params.ke_Central = 3e-2;  % 1/min
end

%% Get Simulated Data
[time, u, y_conc, x_conc] = get_data(true_params, sim_config);

%% Get Fitting Structure
fitObj = FittingStruct(time, u, y_conc, disc_type, elim_type);

%% Run Solver
% Initial Values
V1_init = 1e-2;
V2_init = 1e-3;
CL12_init = 1e-5;
kE_init = 1e-4;

% initial rates
init_vals = struct;
init_vals.kU = 1/V1_init;
init_vals.k12 = CL12_init/V2_init;
init_vals.k21 = CL12_init/V1_init;


init_vals.invsig1 = 1e1;
init_vals.invsig2 = 1e1;

% modify initialization based on elimination type
if strcmp(elim_type, 'constant')
    init_vals.kE = kE_init;
elseif strcmp(elim_type, 'varying')
    init_vals.kE = kE_init*ones(length(fitObj.time)-1, 1);
    init_vals.invsigkE = 1e6;
end

sol_nlp = fitObj.fit_nlp(init_vals);

%% Display Results
% Plot simulated response
fig = figure('Units', 'normalized', 'OuterPosition', [0.1, 0.1, 0.4, 0.8]);
subplot(3,1,1:2)
plot(time, y_conc, '.', 'MarkerSize',8)
hold on
plot(fitObj.time, sol_nlp.x, 'Color', .7*[0,1,0], 'LineWidth',3)
xlim(time([1, end]))

legend('data central','data peripheral','fit')

ylabel('Concentration (\mu M)')
xlabel('Time (min)')

set(gca,'TickDir', 'out', 'box', 'off', 'FontWeight', 'bold')

subplot(3,1,3)
stairs(time, u, 'm', 'LineWidth',2)
ylabel('Injection Rate (\mu mol/min)')

xlim(time([1, end]))
xlabel('Time (min)')

set(gca,'TickDir', 'out', 'box', 'off', 'FontWeight', 'bold')
saveas(fig, 'data_fit', 'svg')

% Show estimated parameters against true parameters if linear elimination
if strcmp(sim_config.model_elim, 'linear')
    T = table;
    T.parameter = fields(true_params);
    T.true = zeros(4, 1);
    T.estimated = zeros(4, 1);
    for i = 1:4
        T.true(i) = true_params.(T.parameter{i});
        T.estimated(i) = sol_nlp.(T.parameter{i});
    end
    disp(T)
end

% Plot estimated elimination rate values
kE = (true_params.Vm_Central./(true_params.Km_Central + x_conc(:,1)))/true_params.Central;

fig = figure;
plot(time, kE,'b', 'LineWidth',2)
hold on
plot(fitObj.time(1:end-1), sol_nlp.ke_Central, 'r-', 'LineWidth',2)
plot(fitObj.time(1:end-1), sol_nlp.ke_Central + sol_nlp.TCI.CI(2:length(fitObj.time)), 'r--')
plot(fitObj.time(1:end-1), sol_nlp.ke_Central- + sol_nlp.TCI.CI(2:length(fitObj.time)), 'r--')
ylabel('k_E (min^{-1})')
xlabel('Time (min)')
legend('True k_E', 'Estimated k_E', '95%% CI')
set(gca,'TickDir', 'out', 'box', 'off', 'FontWeight', 'bold')
saveas(fig, 'kE_fit', 'svg')
