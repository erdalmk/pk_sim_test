classdef FittingStruct
   properties
       time
       inj
       data
       nlp_solver
       type
   end
   methods
      function obj = FittingStruct(t, u, y, type)
          if nargin == 3
              type = 'ForwardEuler';
          elseif nargin~=4
              error('Not enough inputs!')
          end
          % remove duplicate times
          [time, ind_unique, ~] = unique(t);
          obj.time = time;
          obj.inj = u(ind_unique);
          obj.data = y(ind_unique, :);
          obj.nlp_solver = get_solver(length(time), type);
      end
      function sol_struct = fit_nlp(obj, init_struct)
          kU_init = init_struct.kU;
          kE_init = init_struct.kE;
          k12_init = init_struct.k12;
          k21_init = init_struct.k21;
          invsig1_init = init_struct.invsig1;
          invsig2_init = init_struct.invsig2;

          N = length(obj.time);

          % Initial values to pass
          var0 = [kU_init;kE_init;k12_init;k21_init;...
                  invsig1_init;invsig2_init;...
                  zeros(2*(N-1), 1)];
          
          % parameters to pass
          p_val = [obj.time; obj.inj(1:end-1); obj.data(:)];
          
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
          solver = obj.nlp_solver;
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

          sol_struct = struct;
          sol_struct.sol = sol;
          sol_struct.stats = solver.stats;

          sol_struct.kU = kU_sol;
          sol_struct.ke_Central = kE_sol;
          sol_struct.k12 = k12_sol;
          sol_struct.k21 = k21_sol;
          sol_struct.invsig1 = invsig1_sol;
          sol_struct.invsig2 = invsig2_sol;

          sol_struct.Central = 1/kU_sol;
          sol_struct.Peripheral = (k21_sol/kU_sol)/k12_sol;
          sol_struct.Q12 = k21_sol/kU_sol;
          sol_struct.x = x_sol;
      end
   end
end