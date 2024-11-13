%% INIT
clear 
clc
close all
%% SETTINGS
if 1
    %% ENVIRONMENT
    data.env.rho=1025; %water density
    data.env.g=9.80665;%gravity
    data.env.depth=200;%ocean depth
    %% GEOMETRY
    if 1
    %% Locations
    data.geom.locations=[0 0 0].*[1 1 pi/180];% pos_x -- pos_y -- rot_z
    %% Body 1
    data.geom.Body1.D_MC=10;
    data.geom.Body1.D_EC=12.5;
    data.geom.Body1.Draft=20;
    data.geom.Body1.FB=15;
    data.geom.Body1.L_EXT=51.75;
    data.geom.Body1.W_ARM=12.5;
    data.geom.Body1.H_ARM=7;
    data.geom.Body1.W_BR=4;
    data.geom.Body1.H_BR=4;
    data.geom.Body1.t_MC=0.03;
    data.geom.Body1.t_EC=0.0275;
    data.geom.Body1.t_AR=0.0275;
    data.geom.Body1.t_BR=0.025;
    data.geom.Body1.Maxsize_Nemoh=4;
    data.geom.Body1.Minsize_Nemoh=2;
    data.geom.Body1.Platform_Name='VolturnUS';
    data.geom.Body1.PythonScriptName='VolturnUS.py';
    end
    %% TURBINE DATA
    data.WTcomponents.Body1.filepath=[fileparts(pwd) filesep 'mostData' filesep 'windTurbine' ...
        filesep 'turbine_properties' filesep 'Properties_IEA15MW.mat'];

    %% MOORING
    data.mooring.Body1.filepath=[fileparts(pwd) filesep 'mostData' filesep 'mooring' ...
        filesep 'Mooring_IEA15MW_VolturnUS.mat'];
    %% SALOME
    data.salome.SALOME_path=['C:', filesep, 'SALOME-9.9.0'];   
    %% NEMOH
    data.nemoh.frequencies=[68 2*pi*[0.005  1]]; % rad/s
    data.nemoh.wavedir=[1 0 0];
    data.nemoh.symmetric=0; % flag (0/1)
    data.nemoh.results_folder_name='VolturnUS';    
    %% BEMIO
    data.BEMIO.IRF_tEnd=100;
    data.BEMIO.IRF_nDt=1001;
    data.BEMIO.IRF_nDw=1001;
    data.BEMIO.Ex_IRF_tEnd=100;
    data.BEMIO.Ex_IRF_nDt=1001;
    data.BEMIO.Ex_IRF_nDw=1001;
    data.BEMIO.IRF_SS_Omax=6;
    data.BEMIO.IRF_SS_R2t=0.9;
    %% WECSim_source
    data.path.WECSim_source=[fileparts(fileparts(fileparts(pwd))) filesep 'WEC-Sim-MOST'];
    addpath(genpath(data.path.WECSim_source));
    %% FLAGS
    data.flags.plot_hydro_result=1;
    data.flags.plot_mesh_Nemoh=1;
end

%% Main
if 1
    %% Salome
    for i=1:length(fields(data.geom))-1
        %% Write input for salome
        if isfield(data.WTcomponents,['Body' num2str(i)])
            data.WTcomponents.(['Body' num2str(i)])=importdata(data.WTcomponents.(['Body' num2str(i)]).filepath);
            data.geom.(['Body' num2str(i)]).Weight=data.WTcomponents.(['Body' num2str(i)]).m_TOT;
        else
            data.geom.(['Body' num2str(i)]).Weight=0;
        end

        if isfield(data.mooring,['Body' num2str(i)])
            data.mooring.(['Body' num2str(i)]).data=load(data.mooring.(['Body' num2str(i)]).filepath);
            data.geom.(['Body' num2str(i)]).Weight=data.geom.(['Body' num2str(i)]).Weight-...
                data.mooring.(['Body' num2str(i)]).data.F0(3)/data.env.g;
        else
            warning(['No mooring for body ' num2str(i)])
        end

        data.geom.(['Body' num2str(i)]).Rotation_z=data.geom.locations(i,3)*180/pi; %deg
        Platform_data=data.geom.(['Body' num2str(i)]);
        Platform_data_fields=fieldnames(Platform_data);

        fileID = fopen(['SALOME',filesep,'PlatformData.txt'], 'w');
        for ii = 1:length(Platform_data_fields)
            if isnumeric(Platform_data.(Platform_data_fields{ii}))
                fprintf(fileID, '%s\t%f\n', Platform_data_fields{ii}, Platform_data.(Platform_data_fields{ii}));
            elseif ischar(Platform_data.(Platform_data_fields{ii}))
                fprintf(fileID, '%s\t%s\n', Platform_data_fields{ii}, Platform_data.(Platform_data_fields{ii}));
            end
        end
        fclose(fileID);        
        %% Run Salome
        command=['cd "', data.salome.SALOME_path,'" & "',data.salome.SALOME_path,filesep,'W64',filesep,'Python',filesep,'python.exe"' ...
            ' salome -t ' '"' cd filesep 'SALOME' filesep data.geom.(['Body' num2str(i)]).PythonScriptName '"'];
        [stat,~]=system(command);
        if stat~=0
            error('Salome error')
        end
        
        copyfile([cd filesep 'SALOME' filesep 'CAD' filesep data.geom.(['Body' num2str(i)]).Platform_Name '.STEP'],...
            [fileparts(pwd) filesep 'geometry']);
        %% Import data from Salome
        DataOut_Salome=importdata([cd filesep 'SALOME' filesep 'Mass_Inertia' filesep data.geom.(['Body' num2str(i)]).Platform_Name '.dat']);
        rmdir([cd filesep 'SALOME' filesep 'Mass_Inertia'],'s');
        delete([cd filesep 'SALOME' filesep 'PlatformData.txt']);
        clear Platform_data Platform_data_fields
        [Platform.(data.geom.(['Body' num2str(i)]).Platform_Name)]=readSalome(DataOut_Salome,data.geom.locations(i,:));        

        if Platform.(data.geom.(['Body' num2str(i)]).Platform_Name).Negative_BallastMass
            warning('Negative_BallastMass!');
        end
    end
    movefile(['SALOME' filesep 'Mesh'],'NEMOH')
    if not(exist([fileparts(pwd) filesep 'hydroData' ...
            filesep data.nemoh.results_folder_name],'dir'))
        mkdir([fileparts(pwd) filesep 'hydroData' ...
            filesep data.nemoh.results_folder_name]);
    end
    save([fileparts(pwd) filesep 'hydroData' ...
        filesep data.nemoh.results_folder_name filesep 'Mass_Inertia_Properties.mat'],"Platform")
    
    %% Nemoh
    if 1
        %% Pre-proc
        cd('NEMOH')
        assert(FindingNemoh(0,true))        
        %% Run
        NemohRun_from_Salome(Platform,data);
        
        %% h5 file creation
        if 1
        %% Fix Hydrostatic.dat
        pltf_names=fields(Platform);
        for i=1:length(pltf_names)

            lines{1}{1} = sprintf('XF = %f - XG = %f', Platform.(pltf_names{i}).COB(1), Platform.(pltf_names{i}).COG(1));
            lines{1}{2} = sprintf('YF = %f - YG = %f', Platform.(pltf_names{i}).COB(2), Platform.(pltf_names{i}).COG(2));
            lines{1}{3} = sprintf('ZF = %f - ZG = %f', Platform.(pltf_names{i}).COB(3), Platform.(pltf_names{i}).COG(3));
            lines{1}{4} = sprintf('Displacement = %f', Platform.(pltf_names{i}).VolumeDispl);
            lines{1}{5} = sprintf('Waterplane Area = %f', Platform.(pltf_names{i}).A);
            lines{1}{6} = sprintf('Mass = %f', Platform.(pltf_names{i}).mass);

            if length(pltf_names)>1
                fid = fopen([data.nemoh.results_folder_name filesep 'mesh' filesep ...
                    ['Hydrostatics_' num2str(i-1) '.dat']], 'w');
            else
                fid = fopen([data.nemoh.results_folder_name filesep 'mesh' filesep ...
                    'Hydrostatics.dat'], 'w');
            end

            for j = 1:length(lines{1})
                fprintf(fid, '%s\n', lines{1}{j});
            end
            fclose(fid);
        end
        %% KH
        for i=1:length(pltf_names)

            KH=Platform.(pltf_names{i}).Khydro_n_ACOB*data.env.rho*data.env.g;
            KH=KH+[zeros(3,6);
                [zeros(3),Platform.(pltf_names{i}).VolumeDispl*data.env.rho*data.env.g*...
                [-Platform.(pltf_names{i}).COG(3)                 0              Platform.(pltf_names{i}).COG(1);
                0            -Platform.(pltf_names{i}).COG(3)     Platform.(pltf_names{i}).COG(2);
                Platform.(pltf_names{i}).COG(1)            Platform.(pltf_names{i}).COG(2)                0]]];

      
            if length(pltf_names)>1
                fid = fopen([data.nemoh.results_folder_name filesep 'mesh' filesep ...
                    ['KH_' num2str(i-1) '.dat']], 'w');
            else
                fid = fopen([data.nemoh.results_folder_name filesep 'mesh' filesep ...
                    'KH.dat'], 'w');
            end

            for ii = 1:size(KH, 1)
                fprintf(fid, '%E\t', KH(ii,:));
                fprintf(fid, '\n');
            end

            fclose(fid);
        end
        %% Bemio 
        hydro = struct();
        hydro = readNEMOH(hydro,data.nemoh.results_folder_name);
        hydro = radiationIRF(hydro,data.BEMIO.IRF_tEnd,data.BEMIO.IRF_nDt,data.BEMIO.IRF_nDw,[],[]);
        hydro = radiationIRFSS(hydro,data.BEMIO.IRF_SS_Omax,data.BEMIO.IRF_SS_R2t);
        hydro = excitationIRF(hydro,data.BEMIO.Ex_IRF_tEnd,data.BEMIO.Ex_IRF_nDt,data.BEMIO.Ex_IRF_nDw,[],[]);
        writeBEMIOH5(hydro)
        if data.flags.plot_hydro_result
            plotBEMIO(hydro)
        end
        sf=length(findall(0, 'Type', 'figure'));
        spreadfigures(flip(sf:-1:sf-5))

        movefile([data.nemoh.results_folder_name '.h5'],...
            [fileparts(fileparts(pwd)) filesep 'hydroData' filesep data.nemoh.results_folder_name filesep 'hydro.h5'])
        cd ..
        end

    end
end

%% FUNCTIONS
function [Platform]=readSalome(DataOut_Salome,loc)
Platform.location=loc;
% Hydro
headers=string(DataOut_Salome.rowheaders);
Platform.VolumeDispl=DataOut_Salome.data(find(headers=='SubmergedVolume'));
Platform.COB_local=[DataOut_Salome.data(find(headers=='COBx')),DataOut_Salome.data(find(headers=='COBy')),DataOut_Salome.data(find(headers=='COBz'))];
Platform.COB=(Rz(loc(3))*Platform.COB_local')';
Platform.A=DataOut_Salome.data(find(headers=='WaterPlaneArea'));
Platform.Jxx=DataOut_Salome.data(find(headers=='WaterPlaneInertia_Ixx'));
Platform.Jyy=DataOut_Salome.data(find(headers=='WaterPlaneInertia_Iyy'));
Platform.Jxy=DataOut_Salome.data(find(headers=='WaterPlaneInertia_Ixy'));
Platform.Sx=DataOut_Salome.data(find(headers=='WaterPlaneStaticMoment_Sx'));
Platform.Sy=DataOut_Salome.data(find(headers=='WaterPlaneStaticMoment_Sy'));

Platform.Khydro_n_A_local=[...
    0          0           0            0                 0              0 ;
    0          0           0            0                 0              0 ;
    0          0           0       Platform.Sx      -Platform.Sy         0 ;
    0          0           0            0           -Platform.Jxy        0 ;
    0          0           0            0                 0              0 ;
    0          0           0            0                 0              0];
Platform.Khydro_n_A_local=Platform.Khydro_n_A_local+Platform.Khydro_n_A_local'+...
    diag([0    0    Platform.A    Platform.Jxx  Platform.Jyy  0]);

Platform.Khydro_n_A=[Rz(loc(3)) zeros(3);zeros(3) Rz(loc(3))]*Platform.Khydro_n_A_local*[Rz(loc(3))' zeros(3);zeros(3) Rz(loc(3))'];


Platform.Khydro_n_ACOB_local=[...
    0          0           0            0                 0                              0                   ;
    0          0           0            0                 0                              0                   ;
    0          0           0       Platform.Sx      -Platform.Sy                         0                   ;
    0          0           0            0           -Platform.Jxy    -Platform.VolumeDispl*(Platform.COB_local(1)) ;
    0          0           0            0                 0          -Platform.VolumeDispl*(Platform.COB_local(2)) ;
    0          0           0            0                 0                              0                  ];
Platform.Khydro_n_ACOB_local=Platform.Khydro_n_ACOB_local+Platform.Khydro_n_ACOB_local'+...
    diag([0    0    Platform.A    Platform.Jxx+Platform.VolumeDispl*(Platform.COB_local(3))  Platform.Jyy+Platform.VolumeDispl*(Platform.COB_local(3))  0]);

Platform.Khydro_n_ACOB=[Rz(loc(3)) zeros(3);zeros(3) Rz(loc(3))]*Platform.Khydro_n_ACOB_local*[Rz(loc(3))' zeros(3);zeros(3) Rz(loc(3))'];


Platform.Negative_BallastMass=DataOut_Salome.data(find(headers=='Negative_BallastMass'));

% Inertia
groups=["Shell","WaterBallast","SolidBallast"];
Platform.mass=0;
Platform.COG_local=zeros(3,1);
Platform.COG=zeros(3,1);


for i=1:length(groups)

    Platform.bodies.(groups(i)).masses=DataOut_Salome.data(cellfun(@(x) ~isempty(strfind(x, groups(i)) & ...
        strfind(x, '_Mass')), DataOut_Salome.textdata));

    Platform.bodies.(groups(i)).COGs_local=[DataOut_Salome.data(cellfun(@(x) ~isempty(strfind(x, groups(i)) & strfind(x, '_COGx')), DataOut_Salome.textdata))';...
        DataOut_Salome.data(cellfun(@(x) ~isempty(strfind(x, groups(i)) & strfind(x, '_COGy')), DataOut_Salome.textdata))';...
        DataOut_Salome.data(cellfun(@(x) ~isempty(strfind(x, groups(i)) & strfind(x, '_COGz')), DataOut_Salome.textdata))'];
    Platform.bodies.(groups(i)).COGs=Rz(loc(3))* Platform.bodies.(groups(i)).COGs_local;

    Platform.bodies.(groups(i)).mass=sum(Platform.bodies.(groups(i)).masses);

    Platform.bodies.(groups(i)).COG_local=sum(Platform.bodies.(groups(i)).masses'.*Platform.bodies.(groups(i)).COGs_local,2)/Platform.bodies.(groups(i)).mass;
    Platform.bodies.(groups(i)).COG=Rz(loc(3))*Platform.bodies.(groups(i)).COG_local;


    gdl=["_Ixx","_Ixy","_Ixz","_Iyx","_Iyy","_Iyz","_Izx","_Izy","_Izz"];
    for ii=1:3
        for jj=1:3
            kk=(ii-1)*3+jj;
            Platform.bodies.(groups(i)).Inertias_COG_local(ii,jj,:)=reshape(DataOut_Salome.data(cellfun(@(x) ~isempty(strfind(x, groups(i)) & ...
                strfind(x, gdl(kk))), DataOut_Salome.textdata)),1,1,[]);
        end
    end
    for ii=1:size(Platform.bodies.(groups(i)).Inertias_COG_local,3)
        Platform.bodies.(groups(i)).Inertias_COG(:,:,ii)=Rz(loc(3))*Platform.bodies.(groups(i)).Inertias_COG_local(:,:,ii)*Rz(loc(3))';
    end


    Platform.mass=Platform.mass+Platform.bodies.(groups(i)).mass;
    if not(any(isnan(Platform.bodies.(groups(i)).COG)))
        Platform.COG_local=Platform.COG_local+Platform.bodies.(groups(i)).mass*Platform.bodies.(groups(i)).COG_local;
    end


end

Platform.COG_local=Platform.COG_local/Platform.mass;
Platform.COG=Rz(loc(3))*Platform.COG_local;

Platform.I_COG_local=zeros(3,3);
Platform.I_COG=zeros(3,3);

for i=1:length(groups)
    for j=1:size(Platform.bodies.(groups(i)).Inertias_COG,3)

        delta_COG_local=Platform.bodies.(groups(i)).COGs_local(:,j)-Platform.COG_local;
        Platform.I_COG_local=Platform.I_COG_local+...
            Platform.bodies.(groups(i)).Inertias_COG_local(:,:,j)+...
            Platform.bodies.(groups(i)).masses(j)*(delta_COG_local'*delta_COG_local*eye(3)-delta_COG_local*delta_COG_local');
    end
end
Platform.I_SWL_local=Platform.I_COG_local+Platform.mass*(Platform.COG_local'*Platform.COG_local*eye(3)-Platform.COG_local*Platform.COG_local');

Platform.I_COG=Rz(loc(3))*Platform.I_COG_local*Rz(loc(3))';
Platform.I_SWL=Rz(loc(3))*Platform.I_SWL_local*Rz(loc(3))';
end

function [Rz] = Rz(psi)

Rz=[cos(psi)  -sin(psi)  0;
    sin(psi)  cos(psi)   0;
    0          0      1];

end