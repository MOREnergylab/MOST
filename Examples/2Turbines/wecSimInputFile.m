%% Simulation class
simu = simulationClass();                           % Initialize Simulation Class
simu.simMechanicsFile = 'SModel_2Turbines.slx';     % Specify Simulink Model File
simu.mode = 'normal';                               % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer='on';                                 % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                                 % Simulation Start Time [s]
simu.rampTime = 20;                                 % Wave Ramp Time [s]
simu.endTime = 800;                                 % Simulation End Time [s]    
simu.rho = 1025;                                    % Water density [kg/m3]
simu.solver = 'ode4';                               % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.01;                                     % Simulation Time-Step [s]
simu.stateSpace = 0;                                % Flag for convolution integral or state-space calculation, Options: 0 (convolution integral), 1 (state-space)
simu.domainSize = 400;                              % Size of free surface and seabed. This variable is only used for visualization. 
simu.cicEndTime = 60;                               % Specify Convolution integral Time [s]
simu.gravity = 9.80665;                             % Gravity acceleration [m/s2]
simu.b2b = 1;                                       % Flag for body2body interactions, Options: 0 (off), 1 (on)
simu.saveWorkspace=0;                               % Flag to save .mat file for each run, Options: 0 (off), 1 (on)
%% Wave class
% Irregular Waves using Jonswap Spectrum
waves = waveClass('irregular');                     % Initialize WaveClass and Specify Type
waves.phaseSeed = 1;                                % Needed to create different random waves
waves.height = 4;                                   % Significant Wave Height [m]
waves.period = 6;                                   % Peak Period [s]
waves.spectrumType = 'JS';                          % Specify Spectrum Type JS=Jonswap, PM=Pierson-Moskovitz
waves.direction = 0;                                % Wave Directionality [deg]
%% Body class (Platform)
Body_data_folder = fullfile(fileparts(mfilename('fullpath')),'hydroData','2Turbines_Results');
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


body(1).quadDrag.drag = [9.23E+05	0.00E+00	0.00E+00	0.00E+00	-8.92E+06	0.00E+00   % AddBQuad - Additional quadratic drag(N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2) 
                         0.00E+00	9.23E+05	0.00E+00	8.92E+06	0.00E+00	0.00E+00
                         0.00E+00	0.00E+00	2.30E+06	0.00E+00	0.00E+00	0.00E+00
                         0.00E+00	8.92E+06	0.00E+00	1.68E+10	0.00E+00	0.00E+00
                        -8.92E+06	0.00E+00	0.00E+00	0.00E+00	1.68E+10	0.00E+00
                         0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	4.80E+10];

body(2).linearDamping=[100000       0       0       0       0         0   
                            0  100000       0       0       0         0
                            0       0  130000       0       0         0
                            0       0       0       0       0         0
                            0       0       0       0       0         0
                            0       0       0       0       0  13000000];

%% Wind class
wind = windClass();                                                                                     
wind.ConstantWindFlag = 1;                                                                   % Choice of (spatial) constant wind (ConstantWindFlag=1) or (time and spatial) turbolent wind (ConstantWindFlag=0)
wind.V_time_breakpoints=[0  50 600 6000];                                                    % Constant wind speed time breakpoints, the last should be equal or higher than comp. time (only used for ConstantWindFlag=1)
wind.V_modules=[11 11 18 18];                                                                % Constant wind speed modules (only used for ConstantWindFlag=1)
wind.V_directions=[normalize([1,0,0],'norm');...
                   normalize([1,0,0],'norm');...
                   normalize([1,0.1,0],'norm');...
                   normalize([1,0.1,0],'norm')];                                             % Constant wind speed directions (only used for ConstantWindFlag=1)

%% Windturbine class
WT_names={'NREL5MW','IEA15MW'};
aeroLoadsTypes=[0 0];
controls=[0 1];
windSpeed0=11;

for i=1:length(WT_names)

load(fullfile('mostData','windTurbine','control',['SteadyStates_' WT_names{i} '.mat']))
if controls(i)==1; controller='ROSCO';
elseif controls(i)==0; controller='Baseline'; end

windTurbine(i) = windTurbineClass(WT_names{i});                                                                                   % Initialize turbine size and Specify Type
windTurbine(i).aeroLoadsType = aeroLoadsTypes(i);                                                                                 % AeroLoads type: 0-->LUT, 1-->BEM
windTurbine(i).control = controls(i);                                                                                             % Controltype: 0-->Baseline, 1-->ROSCO 
windTurbine(i).omega0 = interp1(SteadyStates.(controller).SS.WINDSPEED,SteadyStates.(controller).SS.ROTSPD,windSpeed0);           % Initial value for rotor speed
windTurbine(i).bladepitch0 = interp1(SteadyStates.(controller).SS.WINDSPEED,SteadyStates.(controller).SS.BLADEPITCH,windSpeed0);  % Initial value for bladepitch
windTurbine(i).GenTorque0= interp1(SteadyStates.(controller).SS.WINDSPEED,SteadyStates.(controller).SS.TORQUE,windSpeed0);        % Initial value for Generator Torque
windTurbine(i).aeroLoadsName = fullfile('mostData','windTurbine','aeroloads',['Aeroloads_' WT_names{i} '.mat']);                  % Aeroloads filename
windTurbine(i).turbineName = fullfile('mostData','windTurbine','turbine_properties',['Properties_' WT_names{i} '.mat']);          % Windturbine properties filename
windTurbine(i).bladeDataName = fullfile('mostData','windTurbine','turbine_properties',['Bladedata_' WT_names{i} '.mat']);         % BladeData filename
windTurbine(i).controlName = fullfile('mostData','windTurbine','control',['Control_' WT_names{i} '.mat']);                        % Controller filename
windTurbine(i).offset_plane=Platform.(pltf_names{i}).location(1:2);                                                               % WindTurbine plane offset with respect w.r.f
windTurbine(i).YawControlFlag = 1;                                                                                                % 0/1 if inactive/active Yaw control
end


%% Constraint
for i=1:length(pltf_names)
constraint(i) = constraintClass(['Constraint' num2str(i)]);              % Initialize constraintClass with name as input
constraint(i).location = [0 0 0];                                        % Constraint Location [m]
end
%% Mooring class
mooring(1) = mooringClass('mooring1');                                                                                            % Initialize mooringClass
mooring(1).location=[Platform.(pltf_names{1}).location(1:2) 0];
mooring(1).lookupTableFile = fullfile(fileparts(mfilename('fullpath')),'mostData','mooring','Mooring_NREL5MW_Spar');              % Load file with mooring look-up table

mooring(2) = mooringClass('mooring2');                                                                     % Initialize mooringClass
mooring(2).location=[Platform.(pltf_names{2}).location(1:2) 0];
mooring(2).Data_moor=struct(...                                                                            % (`string`) Calc Mode ('LUT' or 'NLStatic') Default = ``LUT``
            'd',                                     0.333,...
            'L',                                       850,...        
            'linear_mass_air',                         685,...        
            'number_lines',                              3,...
            'nodes',       [-58   0  -14 ; -837  0  -inf]',...                                             % [Fairlead;Anchor] positions (first line, -inf means water depth)
            'EA',                                   3.27e9,...
            'CB',                                        1,...  
            'MaxIter',                                  150,...
            'TolFun',                                  5e-4,...
            'TolX',                                    5e-4,...
            'HV0_try',                        [1e6    2e6]);
