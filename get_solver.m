function solver = get_solver(N)
    import casadi.*

    % Rate constants
    kU = SX.sym('kU');
    kE = SX.sym('kE');
    k12 = SX.sym('k12');
    k21 = SX.sym('k21');

    % Noise Variances
    invsig1 = SX.sym('invsig1');
    invsig2 = SX.sym('invsig2');

    % Dynamic system parts
    x = SX.sym('x', N-1, 2);
    y = SX.sym('y', N, 2);
    u = SX.sym('u', N-1);
    t = SX.sym('t', N);

    % Construct State Matrix
    A = [-(kE+k21), k21;k12, -k12];
    B = [kU; SX(1,1)];
    
    % Define Dynamics by Forward Euler
    x_all = [SX(1,2); x];
    dxs = A*x_all(1:N-1,:)' + B*u';
    dt = t(2:end)-t(1:end-1);
    x_next = (x_all(1:N-1,:)'+dxs*diag(dt))';
    
    % Pass Dynamics as equality constraints;
    res_eq = x_next-x;
    cons = res_eq(:);

    % Define Log-Likelihood cost
    res = y-x_all;
    J = -N*log(invsig1)/2+invsig1*(res(:,1)'*res(:,1))/2+...
        -N*log(invsig2)/2+invsig2*(res(:,2)'*res(:,2))/2;
    
    % Define variables
    var = [kU;kE;k12;k21;...
           invsig1;invsig2;...
           x(:)];

    % Define parameters
    param = [t;u;y(:)];

    options = struct;
    options.ipopt.max_iter = 1e3;

    nlp = struct('x',var, 'f',J, 'g',cons,'p',param);
    solver = nlpsol('pkfit', 'ipopt', nlp, options);
end