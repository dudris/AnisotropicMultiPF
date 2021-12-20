% input = make_input
% - specifies domain dimensions, initial condition of the phsae fields,
% interface properties of the system including interface width and
% anisotorpy in interface eenrgy and grain boundary mobility
% - values of parameters must be set inside the function
%   - see comments within the function for details
% - 'input' is a structure  containing workspace created in this function
% - from 'input' the phase-field parameters m, kappa, gamma and further time step dt are computed: ' in = input_calc_PFpar_dt(in) '

function input = make_input
%%
Nx = 100; % number of grid points in x direction
Ny = 100; % number of grid points in y direction
precycle = 20; % needed for inclination-dependent simulations. Number of time steps run with isotorpic model to make interface a little diffuse
simtime = 0.01; % simulation time in seconds

IW = 1e-9; % minimal inteface width in meters
IWpts = 7; % number of points in the interface
dx = IW/IWpts; % grid spacing in x direction
dy = dx; % equidistant grid is assumed

% switches
model = 'IWc'; % either of 'IWc', 'IWvG' or 'IWvK', as described in the Paper
is_cond_termd.bool = false; % set true if 

ctrplot = 100; 
ctrcnt = 100; % output made at linspace(1,Ndt,ctrcnt)
plotcond = false;
plotDF = true;
solvermethod = '2Dlin';
laplacianmethod = '9pt20'; % '9pt20', '9pt8', '5pt'
BCs = 'N'; % 'P', 'N', 'mix',

% save phase fields every time after defined number of time steps
outputAtAllCtr = { false , 200}; % should be outputAtAllCtr{2} > ctrplot && mod(outputAtAllCtr{2},ctrplot)==0
% save phase fields after defined simulation time (number of checkpoints between 1st and last time step)
outputAtChckpt = {false, 10}; % number of times in simtime the checkpoint output is saved 
PauseAfterPlotting = false;
PlotAftertstep = Inf; % after selected timestep will plot at every timestep


%% USAGE OF THE FOLLOWING SECTION
% - uncomment the desired initial condition 'ICcond' and set the
% parameters 'ICparam' specifying details of the geometry
% - units used are grid points
% - 2 phases were assumed, denoted as 'solid' and 'liquid'
%   - phase field with number 1 is always 'liquid' 
%   - 'liquid' is represented by a single phase field parameter
% - the following variables must be set 
%   - PFori ... size(PFori) = [1,nOP]
%       - each phase field is defined orientation in polar angle (i.e. w.r.t. x axis)
%       - first phase field assumed to have orientation 0
%   - ind_is_solid ... vector of numbers of phase fields which belong to the 'solid'
%   phase (excluding number 1)
%   - intf ... struct aggregating interface properties
%       - intf.IE_phases ... intetrface energy of distinct types of
%       interfaces in the order [solid-liquid, solid-solid]
%           - systems with liquid-liquid interaction were not simulated,
% PFori = [0,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
% ind_is_solid = [2]; % must not be index 1
% intf.IE_phases = [1/3 ]; % [s-l ]
% intf.is_incl_dep_IE = [ true];
% intf.is_incl_dep_L = [ false ];
%% INITIAL CONDITION 2 PHASE FIELDS
%___ 'CircleInMatrix' ... ICparam = [centerx centery radius]
% ICcode = 'CircleInMatrix';
% ICparam = [Nx/2 , Nx/2 , Nx/3];

% ___ 'Wulff_weak' ... ICparam = [center_x , center_y , radius , Omega], Omega<=1
% ICcode = 'Wulff_weak';
% ICparam = [0.5*Nx 0.5001*Ny Nx/3 0.95];
% ICparam = [0.5*Nx 1.001 Ny*2/3 0.7];

% PFori = [0,0]*(pi/180); % orientation of the PF [1,2,3,...]
% PFphase = [1, 1]; % 
% intf.IE_phases = 1/3 ; 
% intf.mob_phases = 7.5e-16; % m^4/Js
% intf.is_incl_dep_IE = true;
% intf.is_incl_dep_L =  false;

%% INITIAL CONDITION 3 PHASE FIELDS
% ___ 'Tjunction' ... ICparam = [posSS_horizontal , posSL_vertical], 
ICcode = 'Tjunction';
ICparam = [Ny/2 , Nx/2];
% ___ '2CirclesInMatrix' ... ICparam = [center1x center1y radius1 ; center2x center2y radius2]
% ICcode = '2CirclesInMatrix';
% ICparam = [3/4*Nx , 3/4*Nx , Nx/5/sqrt(2) ; Nx/4 Ny/3 Nx/5/sqrt(2)]; 
% ICparam = [0.5*Nx , 0.25*Ny , 0.3*Nx ; 0.5*Nx , 0.75*Ny , 0.3*Nx];

% PFori = [0,10,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
% PFphase = [1, 1, 1]; % 
% intf.IE_phases = [0.3]; % [s-l , s-s], i.e. the types of interfaces 
% intf.is_incl_dep = [false]; % is [s-l , s-s] inclination dependent

PFori = [0,0,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
PFphase = [1, 2, 2]; % 
intf.IE_phases = 0.3*[nan , 1 ; 1 0.4 ]; % [1-1 , 1-2 ; 2-1 , 2-2]
intf.mob_phases = 7.5e-16*ones(2); % m^4/Js
intf.is_incl_dep_IE = false(2);
intf.is_incl_dep_L = false(2);

%% NO USER INPUT (process the interface properties)
nOP = Assign_nOP_from_ICcode(ICcode);
is_inclination_dependent_IE = any(intf.is_incl_dep_IE(:));
is_inclination_dependent_L = any(intf.is_incl_dep_L(:));
is_misori_dependent = length(unique(intf.IE_phases(:)))>1;
%% inclination-dependent properties
if is_inclination_dependent_IE % currently only 1 inclination-dependent type of interface assumed
    intf.params_incl_dep.codeIEaniso = 'IEanisofun_1';
    PFpar_compspec = 'fullaniso'; % 'mean' or 'fullaniso' procesing of IE, IWc is always 'mean'
    intf.isStrongAniso = false; %
    intf.params_incl_dep.nfold = 4 ;
    intf.params_incl_dep.Omega = 0.6; % normalized strength of anisotropy
    intf.params_incl_dep.soaIE = intf.params_incl_dep.Omega/(intf.params_incl_dep.nfold^2-1);
%     intf.params_incl_dep.soaIE = 0.1;
%     intf.params_incl_dep.Omega = intf.params_incl_dep.soaIE*(intf.params_incl_dep.nfold^2-1);
end

if is_inclination_dependent_L % currently only 1 inclination-dependent type of interface assumed
    intf.params_incl_dep.Lnfold = 4 ;
    intf.params_incl_dep.LOmega = 0.2;
    intf.params_incl_dep.LsoaIE = intf.params_incl_dep.LOmega/(intf.params_incl_dep.Lnfold^2-1);
    intf.params_incl_dep.Lfun = @(phi) (ones(size(phi))+intf.params_incl_dep.LsoaIE*cos(phi*intf.params_incl_dep.Lnfold));
    intf.params_incl_dep.override_IW_Lcorr = false; % 'true' to make Lfun the only anisotropy of L
%     intf.params_incl_dep.LsoaIE = 0.1;
%     intf.params_incl_dep.LOmega = intf.params_incl_dep.soaIE*(intf.params_incl_dep.nfold^2-1);
else
    intf.params_incl_dep.override_IW_Lcorr = false; % 
end

% if is_misori_dependent
intf.params_misor_dep.code = 'none' ; % 'none' or 'readshockley'
% intf.params_misor_dep.code = 'readshockley' ; % 'none' or 'readshockley'
% intf.params_misor_dep.width = 15*(pi/180); % rad
% end

%% restrict DF to interace region only
% difference in value of interface-normal magnitude field - the smaller, the wider the interface region
intf.limval.iso = 1e-12; % not expected to be varied much 
% intf_limval_iso = 5e-3; % not expected to be varied much
intf.limval.aniso.outer =  1e-6; %1e-4;
intf.limval.aniso.inner = 5e-3; %1e-2 
% intf.limval.aniso.outer =  1e-4; %;
% intf.limval.aniso.inner = 1e-2; %1e-2 ; ; 
%% check for input consistency

if is_inclination_dependent_IE
    if ~strcmp(PFpar_compspec,'fullaniso')
        warning(['PFpar_compspec is not ''fullaniso''. Expect poorer stability with ''mean''. PFpar_compspec=' PFpar_compspec])
    end
end

validmodels =  categorical({'IWc', 'IWvG','IWvK'});
assert(any(model==validmodels),['Invalid model input: ' model '. Choose one of: ' validmodels])
clear validmodels

% 
numphases = length(unique(PFphase));
assert(all(size(intf.IE_phases)==numphases),['size(intf.IE_phases)~=[numphases,numphases], size(intf.IE_phases)=' num2str(size(intf.IE_phases)) ', numphases=' num2str(numphases)])
assert(all(size(intf.is_incl_dep_IE)==numphases),['size(intf.is_incl_dep_IE)~=[numphases,numphases], size(intf.is_incl_dep_IE)=' num2str(size(intf.is_incl_dep_IE)) ', numphases=' num2str(numphases)])
assert(all(size(intf.is_incl_dep_L)==numphases),['size(intf.is_incl_dep_L)~=[numphases,numphases], size(intf.is_incl_dep_L)=' num2str(size(intf.is_incl_dep_L)) ', numphases=' num2str(numphases)])

%% creating the input structure
save 'dummy.mat'
input = load('dummy.mat');
delete dummy.mat


end%func

%% Assign_nOP_from_ICcode
function nOP = Assign_nOP_from_ICcode(ICcode)
    codes2OP = {'CircleInMatrix','Wulff_weak'};
    codes3OP = {'Tjunction','2CirclesInMatrix'};
    if any(strcmp(ICcode,codes2OP))
        nOP = 2;
    elseif any(strcmp(ICcode,codes3OP))
        nOP =3;
    end
end% func   

