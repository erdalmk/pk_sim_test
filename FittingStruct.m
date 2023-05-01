classdef FittingStruct
   properties
       time
       inj
       data
       nlp_solver
       type
   end
   methods
      function obj = FittingStruct(t, u, y, disc_type, elim_type)
          % remove duplicate times
          [time, ind_unique, ~] = unique(t);
          obj.time = time;
          obj.inj = u(ind_unique);
          obj.data = y(ind_unique, :);
          obj.type = {disc_type, elim_type};
          obj.nlp_solver = get_solver(length(time), disc_type, elim_type);
      end
      function sol_struct = fit_nlp(obj, init_struct)
          % Get Rate initializations
          kU_init = init_struct.kU;
          kE_init = init_struct.kE;
          k12_init = init_struct.k12;
          k21_init = init_struct.k21;
          
          rates_init = [kU_init;kE_init;k12_init;k21_init];
                    
          % Get Variance initializations
          invsig1_init = init_struct.invsig1;
          invsig2_init = init_struct.invsig2;
          sigs_init = [invsig1_init;invsig2_init];
          
          % Parse if kE is varying
          if strcmp(obj.type{2}, 'varying')
            invsigkE_init = init_struct.invsigkE;
          elseif strcmp(obj.type{2}, 'constant')
            invsigkE_init = [];            
          end

          N = length(obj.time);

          % Initial values to pass
          var0 = [rates_init;...
                  sigs_init;...
                  zeros(2*(N-1), 1)];
          
          % parameters to pass
          p_val = [invsigkE_init;obj.time; obj.inj(1:end-1); obj.data(:)];
          
          % constraints to pass
          lbx_val = [1e-1;...
                     1e-6*ones(length(rates_init)-1, 1);...
                     1e-4*ones(length(sigs_init),1);...
                     zeros(2*(N-1), 1)];
          ubx_val = [1e4;...
                     1e0*ones(length(rates_init)-1, 1);...
                     1e3*ones(length(sigs_init),1);...
                     5e2*ones(2*(N-1), 1)];
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
                   
          % get hessian based on laplace integral
          fhess = solver.get_function('nlp_hess_l');
          hess_upper = full(fhess(sol.x, p_val, 1 , sol.lam_g));
          hess_l = hess_upper + triu(hess_upper, 1)';

          fjac_g = solver.get_function('nlp_jac_g');
          [~, jac_g_sparse] = fjac_g(sol.x, p_val);
          jac_gx = full(jac_g_sparse);
          hess_all = [hess_l, jac_gx';...
                      jac_gx, zeros(size(jac_gx, 1))];
          n_params = size(jac_gx, 2) - size(jac_gx, 1);

          mat2sol = [eye(n_params);zeros(size(hess_all, 1)-n_params, n_params)];
          sol_cov = hess_all\mat2sol;
          cov = sol_cov(1:n_params,:);
           
          % save confidence intervals
          sigs = sqrt(diag(cov));
          kEnames = arrayfun(@(x) sprintf('kE%d', x-1), 1:length(kE_init),...
                             'UniformOutput', false)';
          T_sol = table;
          T_sol.names = [{'kU'};kEnames;{'k12';'k21';'invsig1';'invsig2'}];
          T_sol.x = full(sol.x(1:n_params));
          T_sol.CI = 1.96*sigs;
          
          % parse solution
          kU_sol = full(sol.x(1));
          kE_sol = full(sol.x(2:1+length(kE_init)));
          k12_sol = full(sol.x(2+length(kE_init)));
          k21_sol = full(sol.x(3+length(kE_init)));
          invsig1_sol = full(sol.x(4+length(kE_init)));
          invsig2_sol = full(sol.x(5+length(kE_init)));
          x_sol = [zeros(1, 2);...
                   reshape(full(sol.x(6 + +length(kE_init):end)), N-1, 2)];

          sol_struct = struct;
          sol_struct.p = p_val;
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
          
          sol_struct.cov = cov;
          sol_struct.TCI = T_sol;
      end
   end
end