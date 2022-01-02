%% EXAMPLE - run a custom simulation defined in make_input.m
clear 
addpath('input\','solver\')

in = make_input;
% compute model parameters and time step
in = input_calc_PFpar_dt(in);

[p,A,F,in] = run_simulation(in);

%% EXAMPLES - sample simulations from the Paper
% in 'input\examples\' are ready-made input files with sample simulations
% of the four experiments from the Paper
clear
addpath('input\examples\','solver\','input\')

% uncomment the line corresponding to the simulation to be run
% in = input_shrinking_circles; % IWvK, IEtop = 0.3IEbot, IEbot = 0.3 J/m^2
% in = input_trijunction; % IWvG, ratio 1/0.6=1.667
% in = input_wulff; % IWvK, Omega=5, fourfold
% in = input_compensated_aniso; % IWc, Omega=0.6, fourfold


% compute model parameters and time step
in = input_calc_PFpar_dt(in);
% run the simulation
[p,A,F,in] = run_simulation(in);

%% EXAMPLE - only determination of parameters
% see help comments of the functions below for more details
clear 
IEminmax = [0.1, 0.3]; % minimal and maximal interface energy in the system
IEs = [0.3 0.3 0.2 0.1]; % interface energies of the pair-wise interfaces. 
model = 'IWvG'; % either of 'IWc', 'IWvG' or 'IWvK'

GBmobility = 7.5*1e-16*ones(size(IEs)); % interface mobility in m^4/Js
IWmin = 1e-9; % in meters; no interface will be narrower than IWmin

IEinit = determineIEinit(IEminmax,IEs,model);
[kpp0, gam0, m, L, IWout, gsq] = get_PF_parameters(model,IEs, GBmobility ,IWmin, IEinit);
