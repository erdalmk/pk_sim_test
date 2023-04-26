%% Set true parameter values
true_params = struct;
true_params.Central = 5e-2;     % L
true_params.Peripheral = 3e-2;  % L
true_params.Q12 = 3e-4;         % L/min
true_params.ke_Central = 3e-2;  % 1/min

%% Get Simulated Data
[time, u, y_conc, x_conc] = get_data(true_params);

%% Get Nonlinear Solver
N = length(time);
solver = get_solver(N);

%% Run Solver
% Initial Values
V1_init = 1e-2;
V2_init = 1e-3;
CL12_init = 1e-5;

% initial rates
kE_init = 1e-4;
kU_init = 1/V1_init;
k12_init = CL12_init/V2_init;
k21_init = CL12_init/V1_init;

% noise variances
invsig1_init = 1e1;
invsig2_init = 1e1;

% Initial values to pass
var0 = [kU_init;kE_init;k12_init;k21_init;...
        invsig1_init;invsig2_init;...
        zeros(2*(N-1), 1)];

% parameters to pass
p_val = [time; u(1:end-1); y_conc(:)];

% constraints to pass
lbx_val = [1e-1;...
           1e-6*ones(3, 1);...
           1e-3*ones(2,1);...
           zeros(2*(N-1), 1)];
ubx_val = [1e4;...
           1e0*ones(3, 1);...
           1e3*ones(2,1);...
           2.5e2*ones(2*(N-1), 1)];
lbg_val = zeros(2*(N-1), 1);
ubg_val = zeros(2*(N-1), 1);

% run solver
sol = solver('x0',var0,...
             'p',p_val,...
             'lbx',lbx_val,...
             'ubx',ubx_val,...
             'lbg',lbg_val,...
             'ubg',ubg_val);

% parse solution
kU_sol = full(sol.x(1));
kE_sol = full(sol.x(2));
k12_sol = full(sol.x(3));
k21_sol = full(sol.x(4));
invsig1_sol = full(sol.x(5));
invsig2_sol = full(sol.x(6));
x_sol = [zeros(1, 2);...
         reshape(full(sol.x(7:end)), N-1, 2)];

%% Plot simulated response
figure
plot(time, y_conc, '.')
hold on
set(gca,'ColorOrderIndex',1)
plot(time, x_conc)
plot(time, x_sol)
ylabel('Concentration (\mu M)')

yyaxis right
stairs(time, u, 'r')
ylabel('Injection Rate (\mu mol/min)')

xlim(time([1, end]))
xlabel('Time (min)')

ax = gca;
ax.YAxis(2).Color = 'r';
set(ax,'TickDir', 'out', 'box', 'off')