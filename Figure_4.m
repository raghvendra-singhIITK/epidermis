% Initial conditions
% [Notch_TAC, Wnt_SC, Wnt_TAC, YAP_TAC]
X0 = [0; 0; 0; 0];  

% Rate constants
k = [1; 0.2; 0.05;      % Notch  (k_Notch, kd_Notch, kp_Notch)
     1.4; 0.2; 0.1;     % Wnt_SC (k_WntSC, kd_WntSC, kp_WntSC)
     1; 0.2; 0.1;       % Wnt_TA (k_WntTA, kd_WntTA, kp_WntTA)
     0.188; 0.02; 0.1]; % YAP    (k_YAP,   k_YAP,    k_YAP)

% Hill coefficient
n = 2;

% Run time
Tend = 1000;
% Step size
dt = 0.1;
% Time grid
t_uniform = 0:dt:Tend;
num_timepoints = length(t_uniform);

% Number of simulations
numSims = 1000;

% Initial perturbation
perturbation_size = 1e-2;
% Continious Noise
noise_strength    = 1e-3;  
% Fix random seed for reproducing results
rng(1);

% PARTIAL STOCHASTIC SIMULATIONS BASED ON EULER-MARUYAMA METHOD

% Stores all simulations
all_X = zeros(4, num_timepoints, numSims);

% Adds small random perturbations to initial conditions
for sim = 1:numSims
    delta = perturbation_size * randn(4,1);
    X = max(X0 + delta, 0);
    X_traj = zeros(4, num_timepoints);
    X_traj(:,1) = X;

    for tidx = 2:num_timepoints
        t = t_uniform(tidx);
        f = odes_epidermal(X, k, n);
        % Brownian increment
        dW = sqrt(dt) * randn(4,1);  
        X = X + f*dt + noise_strength * dW;
        X = max(X,0); 
        X_traj(:,tidx) = X;
    end
    all_X(:,:,sim) = X_traj;
end

% Mean and Standard deviation of Stochastic results
mean_X = mean(all_X, 3);
std_X  = std(all_X, 0, 3);

% Reference unperturbed dynamic system
[~, X_ode] = ode45(@(t,x) odes_epidermal(x,k,n), t_uniform, X0);
X_ode = X_ode';  

% Plots
species_names = {'Notch activity in the TAC, A_{Notch }', ...
                 'Wnt activity in the SC, A_{WntSC} ', ...
                 'Wnt activity in the TAC, A_{WntTAC} ', ...
                 'YAP activity in the TAC, A_{YAP} '};
colors = lines(4);

for s = 1:4
    figure;
    hold on;
    % plot shaded region for mean ± standard deviation
    fill([t_uniform fliplr(t_uniform)], ...
         [mean_X(s,:) + std_X(s,:), fliplr(mean_X(s,:) - std_X(s,:))], ...
         colors(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % plot mean of stochastic
    plot(t_uniform, mean_X(s,:), '-', 'Color', colors(s,:), 'LineWidth', 4);
    % plot deterministic reference
    plot(t_uniform, X_ode(s,:), '--', 'Color', [0 0 0], 'LineWidth', 3);
    
    % Axis labels and title
    xlabel('Time (a.u.)', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel([species_names{s} ' (a.u.)'], 'FontName', 'Times New Roman', 'FontSize', 20);
       
    % Legend
    lgd = legend('The Stochastic trajectories','Average of the Stochastic trajectories','Deterministic Dynamics','Location','Best');
    set(lgd, 'FontName', 'Times New Roman', 'FontSize', 20, 'Box','off');
    
    % Grid and styling
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 26, 'LineWidth', 1.5);
 end

% Functions

% Time derivative of state vector
function dxdt = odes_epidermal(x, k, n)
    % % Reaction rates (propensities)
    a = propensities(x, k, n);
    dxdt = zeros(4,1);
    % ODEs describing Regulation, Degradation, Generation terms
    dxdt(1) = a(1) - a(2) + a(3);    % Notch_TAC
    dxdt(2) = a(4) - a(5) + a(6);    % Wnt_SC
    dxdt(3) = a(7) - a(8) + a(9);    % Wnt_TAC
    dxdt(4) = a(10) - a(11) + a(12); % YAP_TAC
end

% Calculates propensities
function a = propensities(X, k, n)
    Notch_TAC = X(1);
    Wnt_SC    = X(2);
    Wnt_TAC   = X(3);
    YAP_TAC   = X(4);

    WSCn  = Wnt_SC^n;
    WTACn = Wnt_TAC^n;
    YTACn = YAP_TAC^n;
    NTACn = Notch_TAC^n;

    a = zeros(12,1);
    % Notch_TAC propensities
    a(1)  = k(1)  * (WSCn / (1 + WSCn + WTACn + YTACn)); 
    a(2)  = k(2)  * Notch_TAC;
    a(3)  = k(3);
    % Wnt_SC propensities
    a(4)  = k(4)  * (NTACn / (1 + NTACn));
    a(5)  = k(5)  * Wnt_SC;
    a(6)  = k(6);
    % Wnt_TAC propensities
    a(7)  = k(7)  * (1 / (1 + NTACn));
    a(8)  = k(8)  * Wnt_TAC;
    a(9)  = k(9);
    % YAP_TAC propensities
    a(10) = k(10) * (NTACn / (1 + NTACn));
    a(11) = k(11) * YAP_TAC;
    a(12) = k(12);
end
