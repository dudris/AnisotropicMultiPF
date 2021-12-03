% input = inputfile_GGaniso
%   - to load parameters into Multip_aniso 
%   - values of parameters must be set inside the function
%   - 'input' is a struct containing workspace created in this function
%   - from 'input' PF parameters and dt are computed: ' in = input_calc_PFpar_dt(varargin) '

function input = inputfile_GGaniso
%%
Nx = 100;
Ny = 100;
precycle = 100;
simtime = 0.03;
% simtime = 1e5;
Ndtstart = 1;

% IW = 1e-9/(Nx/100);
IW = 1e-9;
% IW = 1e-6;
dx = IW/7;
dy = dx;

% switches
model = 'IWvG'; % IWc, IWvG or IWvK
% is_with_constant_IW = false;
is_conc_conserved = false;
is_conserved = false;
is_cond_termd.bool = false; 
spec_singlePF_IEcalc = 'default'; %  either 'default' or PF number to calculate IE of. Option 'default' sets the nOP-th grain
is_fixedPF.bool = false; % 

ctrplot = 100; 
ctrcnt = 100; % output made at linspace(1,Ndt,ctrcnt)
plotcond = true;
plotDF = false;
solvermethod = '2Dlin';
laplacianmethod = '9pt20'; % '9pt20', '9pt8', '5pt'
BCs = 'N'; % 'P', 'N', 'mix', 'Nspecial'
% BCs = 'Nspecial'; % 'P', 'N', 'mix', 'Nspecial'
outputAtAllCtr = { false , 200}; % should be outputAtAllCtr{2} > ctrplot && mod(outputAtAllCtr{2},ctrplot)==0
outputAtChckpt = {false, 10}; % number of times in simtime the checkpoint output is saved 
PauseAfterPlotting = false;
PlotAftertstep = Inf; % after selected timestep will plot at every timestep
save_OPs_after_sim = false; % to save all output PFs to results file
Ndtstart_is_1 = true;
% unused switches 
all_expr_anal = false;
calc_arc_length = false;
dispDFdetailscond = false;

if is_cond_termd.bool
    is_cond_termd.code = 'area_ratio';
    is_cond_termd.PFnum = 2;
    is_cond_termd.ratio = 0.6;
%     is_cond_termd.ratio = (1/3)/0.38;
    
%     is_cond_termd.code = 'area_abs';
%     is_cond_termd.PFnum = 2;
%     is_cond_termd.area = 5e-17;
%     is_cond_termd.area = pi*(2.5e-9)^2/2; % in m^2, area of semidisc of radius 2.5 nm
    
%     is_cond_termd.code = '1g_meanIE'; % '1g_totIE' or '1g_meanIE'
%     is_cond_termd.PFnum = spec_singlePF_IEcalc;
% %     is_cond_termd.rellim = 1e-6;
% %     is_cond_termd.rellim = 5e-5; % sufficient for 1g_totIE
%     is_cond_termd.rellim = 5e-5; % sufficient for 1g_meanIE
%     is_cond_termd.meantime = simtime/50; % s , time to take mean of energy in for assessing term. cond.
    
end

if is_fixedPF.bool
    is_fixedPF.code = 'mobility';
    is_fixedPF.PFnum = 2;
    is_fixedPF.factor = 0.5e-2;
end

if strcmp(BCs,'Nspecial')
% right-half-plane angle must be with -sign
%         BCs_specs.bottom = [1, 50, 120 ; 51, Nx, -30]; % [ind11, ind12, angle1 ; ind21, ind22, angle ; ind31 ...]
%     BCs_specs.top = [1, 50, 120 ; 51, Nx, -30]; % [ind11, ind12, angle1 ; ind21, ind22, angle ; ind31 ...]

%     BCs_specs.bottom = [1, Nx, 175 ]; % [ind11, ind12, angle1 ; ind21, ind22, angle ; ind31 ...]
%     BCs_specs.bottom = [1, Nx, NN(1)]; % [ind11, ind12, angle1 ; ind21, ind22, angle ; ind31 ...]
%     BCs_specs.top = [1, Nx, -BCs_specs.bottom(3)]; % tiled plane simulation
    
    BCs_specs.bottom = [1, floor(Nx/2), 120 ; ceil(Nx/2), Nx, -100]; % [ind11, ind12, angle1 ; ind21, ind22, angle ; ind31 ...]
    BCs_specs.top = [1, floor(Nx/2), 180 ; ceil(Nx/2), Nx, -180];
    BCs_specs.left = [1, floor(Ny/2), 60 ; ceil(Ny/2), Ny, -45];
    BCs_specs.right = [1, floor(Ny/2), 45 ; ceil(Ny/2), Ny, -60];
end

%% model input assertion
validmodels =  categorical({'IWc', 'IWvG','IWvK'});
assert(any(model==validmodels),['Invalid model input: ' model '. Choose one of: ' validmodels])
clear validmodels
%% IC 2 PFs 
%___ 'CircleInMatrix' ... ICparam = [centerx centery radius]
% ICcode = 'CircleInMatrix';
% ICparam = [Nx/2 , 0.3*Ny , 0.38*Nx];
% ICparam = [Nx/2 , 0.3*Ny , Nx/3];
% ICparam = [Nx/2 , Nx/2 , Nx/3];
% ICparam = [0 , 0 , Nx/2];
% ICparam = [Nx/2 , 0 , Nx/3];
% ICparam = [Nx/2, 1 , Ny/3*2];
%___ 'VerticalPlane' ... ICparam = [ position ] , e.g. [Nx/2]
% ICcode = 'VerticalPlane';
% ICparam = Nx/2;
%___ 'TiltedPlane' ... ICparam = [ y1intercept yNyintercept ]  ... straight line intersecting y=1 and y=Ny in points param(1) and param(2)
% ICcode = 'TiltedPlane';
%  ICparam = [Nx/4 3*Nx/4];
%___ 'TiltedPlaneAng' ... ICparam = [ inclination_angle_degrees ]  ... straight line passing through [Nx/2, Ny/2] tilted under 'inclination_angle_degrees'
% ICcode = 'TiltedPlaneAng';
% ICparam = 60;
%___ 'VerticalSlab' ... ICparam = [position, halfwidth]
% ICcode = 'VerticalSlab';
%  ICparam = [ceil(Nx/2) ceil(Nx/7)];
%___ 'EllipseInMatrix' ... ICparam = [centerx , centery , semiaxX, semiaxY]
% ICcode = 'EllipseInMatrix';
%  ICparam = [Nx/2 , Nx/2 , IW/3 , IW/10];
%  ICparam = [Nx/2 , 1 , Nx/3 , Ny/2];
%___ 'RectangleInMatrix' ... ICparam = [half_edge_x half_edge_y]
% ICcode = 'RectangleInMatrix';
% ICparam = [Nx/4 , Ny/4];
% ___ 'Wulff_weak' ... ICparam = [center_x , center_y , radius , Omega]
ICcode = 'Wulff_weak';
ICparam = [0.5*Nx 0.5001*Ny Nx/3 0.95];
% ICparam = [0.5*Nx 1.001 Ny*2/3 0.7];
%___ 'EAVarAng' ... ICparam = [angLdeg angRdeg init_rad_meters ]
    % compound of 2 circular arcs with different contact angle with x axis connected by tangent  line
% ICcode = 'EAVarAng'; % Energy Assessment with Variable Angle                        
% ICparam = [45 60, 3.33779e-09];
% ICparam = [45 60, 50*dx];

PFori = [0,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
ind_is_solid = [2]; % must not be index 1
intf.IE_phases = [1/3 ]; % [s-l ]
intf.is_incl_dep_IE = [ true];
intf.is_incl_dep_L = [ false ];

%% IC 3 PFs 
%___ '3junctions' ... ICparam = [3jun_pos , angle], 3jun_pos<=Ny/2 
% ICcode = '3junctions';
% ICparam = [Ny/3 , pi/4]; 
% ___ 'SemiCircleOnPlane' ... ICparam = [plane_pos , radius], 
% ICcode = 'SemiCircleOnPlane';
% ICparam = [0.1*Ny , Nx/3];
% ___ 'Tjunction' ... ICparam = [posSS_horizontal , posSL_vertical], 
% ICcode = 'Tjunction';
% ICparam = [Ny/2 , Nx/2];
% ___ '2CirclesInMatrix' ... ICparam = [center1x center1y radius1 ; center2x center2y radius2]
% ICcode = '2CirclesInMatrix';
% ICparam = [3/4*Nx , 3/4*Nx , Nx/5/sqrt(2) ; Nx/4 Ny/3 Nx/5/sqrt(2)]; 
% ICparam = [0.5*Nx , 0.25*Ny , 0.3*Nx ; 0.5*Nx , 0.75*Ny , 0.3*Nx];

% PFori = [0,10,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
% ind_is_solid = [1,2,3]; % PF indices which are 'solid'
% intf.IE_phases = [0.3]; % [s-l , s-s], i.e. the types of interfaces 
% intf.is_incl_dep = [false]; % is [s-l , s-s] inclination dependent
% 

% PFori = [0,0,0]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
% ind_is_solid = [2,3]; % PF indices which are 'solid'
% intf.is_incl_dep_IE = [false false]; % is [s-l , s-s] inclination dependent 
% intf.is_incl_dep_L = [ true false];
% intf.IE_phases = 0.3*[1 , 0.4]; % [s-l , s-s], i.e. the types of interfaces 
% intf.is_incl_dep_IE = [false false]; % is [s-l , s-s] inclination dependent
% intf.is_incl_dep_L = [ false false];
%% IC 4 PFs 
% ___ 'SemiCircleOnTjunction' ... ICparam = [plane_pos , radius], 
% ICcode = 'SemiCircleOnTjunction';
% ICparam = [Ny/3 , Nx/3]; 
% ___ 'TestEnergyCalc_4PFs' ... ICparam = nan, Nx >= 15*IW/dx, Ny >=9*IW/dx
% ICcode = 'TestEnergyCalc_4PFs';
% ICparam = nan; % uncomment for run
% ___'3CirclesInMatrix' ... ICparam = [center1x center1y radius1 ; center2x center2y radius2 ; center3x center3y radius3]
% ICcode = '3CirclesInMatrix';
% ICparam = [0.3*Nx , 0.25*Ny , 0.2*Nx ; 0.3*Nx , 0.75*Ny , 0.2*Nx, ; 0.7*Nx , 0.5*Ny , 0.2*Nx];

% 
% PFori = [0,10,20,30]*(pi/180); % orientation of the PF [1,2,3,...], liquid has always 0
% ind_is_solid = [2,3,4]; % PF indices which are 'solid'
% intf.IE_phases = 0.3*[1 , 0.6]; % [s-l , s-s], i.e. the types of interfaces 
% intf.is_incl_dep_IE = [false false]; % is [s-l , s-s] inclination dependent 
% intf.is_incl_dep_L = [ false false];
%% process the interface properties 
nOP = Assign_nOP_from_ICcode(ICcode);
is_inclination_dependent_IE = any(intf.is_incl_dep_IE);
is_inclination_dependent_L = any(intf.is_incl_dep_L);
is_misori_dependent = length(intf.IE_phases)>1;

% if is_misori_dependent
intf.params_misor_dep.code = 'none' ; % 'none' or 'readshockley'
% intf.params_misor_dep.code = 'readshockley' ; % 'none' or 'readshockley'
% intf.params_misor_dep.width = 15*(pi/180); % rad
% end

if is_inclination_dependent_IE % currently only 1 inclination-dependent type of interface assumed
    PFpar_compspec = 'fullaniso'; % 'mean' or 'fullaniso' procesing of IE, IWc is always 'mean'
    intf.isStrongAniso = false;
    intf.params_incl_dep.nfold = 4 ;
    intf.params_incl_dep.Omega = 3.6;
    intf.params_incl_dep.soaIE = intf.params_incl_dep.Omega/(intf.params_incl_dep.nfold^2-1);
%     intf.params_incl_dep.soaIE = 0.1;
%     intf.params_incl_dep.Omega = intf.params_incl_dep.soaIE*(intf.params_incl_dep.nfold^2-1);
    intf.params_incl_dep.codeIEaniso = 'IEanisofun_1';
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

% % only 1 type of intetrface assumed
% isLaniso = false;
% if isLaniso
%     delta = 0.05;
%     Lfun = @(phi) (ones(size(phi))+delta*cos(phi*nfold));
% end


%% concentration volume conservation specification
if is_conc_conserved
    PF_to_conserve = { [2] , []}; % liquid phase primarily conserved (PFnum 1), 2nd PF in 'PF_to_conserve' can be any other
        assert(all(ismember([PF_to_conserve{:}],1:nOP)),['is_conc_conserved=true and PF_to_conserve=[' sprintf('%1.0f,', [PF_to_conserve{:}]) '] but nOP=' num2str(nOP)])
    conserve_2PFs = ~isempty(PF_to_conserve{2});
    pref_GP_diff = AssignGPdifferencePrefactor(nOP,PF_to_conserve);
    cLeq = 0.98;%0.98; % equilibrium molar fraction in liquid
    cSeq = 0.02 ; %0.02; % equilibrium molar fraction in solid
    AL = 1*10^10.5; %[J/m3] prefactor for parabolic fit of free energy
%     D = 10^(-15.35);
%     D =  0.3*(dx^2)/dt;
    D =  0.3*(dx^2)/dt;
    Mconc_L = D/AL; % Vishal : 10^(-16)/AL ~ 10^(-23)
    % with stabcritC = D*dt/(dx^2) < 0.18 it is stable
%     D*dt/(dx^2);
%     Mconc_L = 10^(-27); % Vishal : 10^(-16)/AL ~ 10^(-23)
end

%% restrict DF to interace region only
% difference in value of interface-normal magnitude field - the smaller, the wider the interface region
intf.limval.iso = 1e-12; % not expected to be varied much 
% intf_limval_iso = 5e-3; % not expected to be varied much
intf.limval.aniso.outer =  1e-6; %1e-4;
intf.limval.aniso.inner = 5e-3; %1e-2 
% intf.limval.aniso.outer =  1e-4; %;
% intf.limval.aniso.inner = 1e-2; %1e-2 ; ; 
%% check for input consistency
single_vol_cons_approach = (is_conc_conserved & ~is_conserved) | (~is_conc_conserved & is_conserved) |  (~is_conc_conserved & ~is_conserved);
% true when only one of the switches is true and the other false.
assert(single_vol_cons_approach,'inputfile_GGaniso: Modify input. Volume conservation by both Lagr. multipliers and concentration field applied.')
clear single_vol_cons_approach ans

if nOP == 2
    intf_props = length(ind_is_solid) == 1;
    assert(intf_props,'For nOP==2 There should be only 1 interface but there are probably more specified.')
    clear intf_props
end 

if ischar(spec_singlePF_IEcalc) % string input
    assert(strcmp(spec_singlePF_IEcalc,'default'),['Invalid specifyer spec_singlePF_IEcalc. Only ''default'' accepted when class char.'])
elseif isnumeric(spec_singlePF_IEcalc)
    assert(any((1:nOP)-spec_singlePF_IEcalc==0),'Invalid specifyer spec_singlePF_IEcalc. Only spec_singlePF_IEcalc>=1 & spec_singlePF_IEcalc<=nOP accepted when class numeric.')
end

if strcmp(BCs,'Nspecial')
    assert(nOP ==2,'BCs=Nspecial but nOP~=2.')
    if ~strcmp(model,'IWvK')
        warning(['Inconsistent input: model= ' model ' and BCs=Nspecial. BCs changed to N.'])
        BCs = 'N';
    end
end

if is_inclination_dependent_IE
    if ~strcmp(PFpar_compspec,'fullaniso')
        warning(['PFpar_compspec is not ''fullaniso''. Expect poorer stability with ''mean''. PFpar_compspec=' PFpar_compspec])
    end
end

%% creating the input structure
save 'dummy.mat'
input = load('dummy.mat');
delete dummy.mat


end%func

%% Assign_nOP_from_ICcode
function nOP = Assign_nOP_from_ICcode(ICcode)
    codes2OP = {'CircleInMatrix','VerticalPlane','LoadContour','TiltedPlane','TiltedPlaneAng','VerticalSlab','EllipseInMatrix','RectangleInMatrix','Wulff_weak','EAVarAng'};
    codes3OP = {'3junctions','SemiCircleOnPlane', 'SemiCircleOnPlaneDiffuse' ,'Tjunction','2CirclesInMatrix'};
    codes4OP = {'SemiCircleOnTjunction','SemiCircleOnTjunctionDiffuse','TestEnergyCalc_4PFs','3CirclesInMatrix'};
    if any(strcmp(ICcode,codes2OP))
        nOP = 2;
    elseif any(strcmp(ICcode,codes3OP))
        nOP =3;
    elseif any(strcmp(ICcode,codes4OP))
        nOP =4;
    end
end% func   

%% AssignGPdifferencePrefactor
function pref_GP_diff = AssignGPdifferencePrefactor(nOP,PF_to_conserve)
    pref_GP_diff = -ones(nOP,2);
    pref_GP_diff(PF_to_conserve{1},1) =1;
    pref_GP_diff(PF_to_conserve{2},2) =1;
end

