%% Simulation class
simu = simulationClass();                           % Initialize Simulation Class
simu.simMechanicsFile = 'SModel_HybridSpar.slx';    % Specify Simulink Model File
simu.mode = 'normal';                               % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer='on';                                 % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                                 % Simulation Start Time [s]
simu.rampTime = 20;                                 % Wave Ramp Time [s]
simu.endTime = 800;                                 % Simulation End Time [s]    
simu.rho = 1025;                                    % Water density [kg/m3]
simu.solver = 'ode4';                               % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.01;                                     % Simulation Time-Step [s]
simu.stateSpace = 0;                                % Flag for convolution integral or state-space calculation, Options: 0 (convolution integral), 1 (state-space)
simu.domainSize = 100;                              % Size of free surface and seabed. This variable is only used for visualization. 
simu.cicEndTime = 60;                               % Specify Convolution integral Time [s]
simu.gravity = 9.80665;                             % Gravity acceleration [m/s2]
simu.b2b = 1;                                       % Flag for body2body interactions, Options: 0 (off), 1 (on). Default = ``0``
simu.saveWorkspace=0;                               % Flag to save .mat file for each run, Options: 0 (off), 1 (on). Default = ``1``
%% Wave class
% Irregular Waves using Jonswap Spectrum
waves = waveClass('irregular');                     % Initialize WaveClass and Specify Type
waves.phaseSeed = 1;                                % Needed to create different random waves
waves.height = 1;                                   % Significant Wave Height [m]
waves.period = 6;                                   % Peak Period [s]
waves.spectrumType = 'JS';                          % Specify Spectrum Type JS=Jonswap, PM=Pierson-Moskovitz
waves.direction = 0;                                % Wave Directionality [deg]
%% Body class (Platform)
Body_data_folder = fullfile(fileparts(mfilename('fullpath')),'hydroData','HybridSpar');
load([Body_data_folder filesep 'Mass_Inertia_Properties.mat'])
pltf_names=fields(Platform);

for i=1:length(pltf_names)
body(i) = bodyClass([Body_data_folder filesep 'hydro.h5']);                                    %#ok<*SAGROW> % Initialize bodyClass (giving hydro data file as input)
body(i).geometryFile = ['geometry' filesep  pltf_names{i} '.STEP'];                            % Geometry File 
body(i).mass = Platform.(pltf_names{i}).mass;                                                  % User-Defined mass [kg]
body(i).inertia = diag(Platform.(pltf_names{i}).I_COG);                                        % Moment of Inertia (diagonal part) [kg-m^2]
body(i).inertiaProducts=Platform.(pltf_names{i}).I_COG([4 7 8]);                               % Moment of Inertia (extradiagonal part) [kg-m^2]
body(i).initial.displacement=[Platform.(pltf_names{i}).location(1:2) 0];
end


body(1).linearDamping=[150000       0       0       0       0         0   
                            0  150000       0       0       0         0
                            0       0  200000       0       0         0
                            0       0       0       0       0         0
                            0       0       0       0       0         0
                            0       0       0       0       0  13000000];

%% Wind class
wind = windClass();                                                                                     
wind.ConstantWindFlag = 1;                                                             % Choice of (spatial) constant wind (ConstantWindFlag=1) or (time and spatial) turbolent wind (ConstantWindFlag=0)
wind.V_time_breakpoints=[0 6000];                                                      % Constant wind speed time breakpoints, the last should be equal or higher than comp. time (only used for ConstantWindFlag=1)
wind.V_modules=[18 18];                                                                % Constant wind speed modules (only used for ConstantWindFlag=1)
wind.V_directions=[normalize([1,0,0],'norm');...
                   normalize([1,0,0],'norm')];                                         % Constant wind speed directions (only used for ConstantWindFlag=1)

%% Windturbine class
windSpeed0=18;
load(fullfile('mostData','windTurbine','control','SteadyStates_IEA15MW.mat'))

windTurbine(1) = windTurbineClass('IEA15MW');                                                                                         % Initialize turbine size and Specify Type
windTurbine(1).aeroLoadsType = 0;                                                                                                     % AeroLoads type: 0-->LUT, 1-->BEM
windTurbine(1).control = 1;                                                                                                           % Controltype: 0-->Baseline, 1-->ROSCO 
windTurbine(1).omega0 = interp1(SteadyStates.ROSCO.SS.WINDSPEED,SteadyStates.ROSCO.SS.ROTSPD,windSpeed0,"linear","extrap");           % Initial value for rotor speed
windTurbine(1).bladepitch0 = interp1(SteadyStates.ROSCO.SS.WINDSPEED,SteadyStates.ROSCO.SS.BLADEPITCH,windSpeed0,"linear","extrap");  % Initial value for bladepitch
windTurbine(1).GenTorque0= interp1(SteadyStates.ROSCO.SS.WINDSPEED,SteadyStates.ROSCO.SS.TORQUE,windSpeed0,"linear","extrap");        % Initial value for Generator Torque
windTurbine(1).aeroLoadsName = fullfile('mostData','windTurbine','aeroloads','Aeroloads_IEA15MW.mat');                                % Aeroloads filename
windTurbine(1).turbineName = fullfile('mostData','windTurbine','turbine_properties','Properties_IEA15MW.mat');                        % Windturbine properties filename
windTurbine(1).bladeDataName = fullfile('mostData','windTurbine','turbine_properties','Bladedata_IEA15MW.mat');                       % BladeData filename
windTurbine(1).controlName = fullfile('mostData','windTurbine','control','Control_IEA15MW.mat');                                      % Controller filename
windTurbine(1).offset_plane=Platform.(pltf_names{1}).location(1:2);                                                                   % WindTurbine plane offset with respect w.r.f

%% Constraint
constraint(1) = constraintClass('Constraint1');              % Initialize constraintClass with name as input
constraint(1).location = [0 0 0];                            % Constraint Location [m]

%% Mooring class
mooring(1) = mooringClass('mooring1');                                                                                            % Initialize mooringClass
mooring(1).location=[Platform.(pltf_names{1}).location(1:2) 0];
mooring(1).lookupTableFile = fullfile(fileparts(mfilename('fullpath')),'mostData','mooring','Mooring_NREL15MW_Spar');              % Load file with mooring look-up table

%% PTO Class
% Translational PTOs
pto(1) = ptoClass('PTO1');              	
pto(1).stiffness=0;                             	
pto(1).damping=1200000;                       	
pto(1).location = [Platform.(pltf_names{i}).location(1:2) 0]+Platform.(pltf_names{2}).COG';   
