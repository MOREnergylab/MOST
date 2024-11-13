%--------------------------------------------------------------------------------------
%
% --> function found = FindingNemoh(QTF, verbose)
%
% Purpose: Check that the Nemoh executables are in the system path, add
% them otherwise.
%
% Inputs :
% - QTF    : Whether to also include the QTF programs in the search
% - verbose: Whether to output information about the search in the console
%

function found = FindingNemoh(QTF, verbose)

if nargin<1 
    QTF = true;
end
if nargin<2
    verbose = false;
end

found = false;

if isunix
    search_command = 'which';
    executable_suffix = '';
else
    search_command = 'where';
    executable_suffix = '.exe';
end

% Search system PATH
[mesh_status, ~] = system([search_command ' mesh']);
[hydrosCal_status, ~] = system([search_command ' hydrosCal']);
[preProc_status, ~] = system([search_command ' preProc']);
[solver_status, ~] = system([search_command ' solver']);
[postProc_status, ~] = system([search_command ' postProc']);
[QTFpreproc_status, ~] = system([search_command ' QTFpreProc']);
[QTFsolver_status, ~] = system([search_command ' QTFsolver']);
[QTFpostProc_status, ~] = system([search_command ' QTFpostProc']);
status = all([mesh_status, hydrosCal_status, preProc_status, solver_status, postProc_status]);
status_QTF = all([QTFpreproc_status, QTFsolver_status, QTFpostProc_status]);
if status
    if verbose
        disp('Nemoh programs not found in system PATH');
    end
elseif ~QTF % Basic Nemoh found and QTF not needed
    if verbose
        disp('Nemoh programs found in system PATH');
    end
    found = true;
    return
elseif status_QTF % QTF programs needed but not found
    if verbose
        disp('Nemoh programs found in system PATH but not for QTF');
    end
else % Basic Nemoh and QTF found
    if verbose
        disp('Nemoh programs (including QTF) found in system PATH');
    end
    found = true;
    return
end

% Search likely locations
[pathstr,~,~] = fileparts(mfilename('fullpath'));
nemoh_root = [pathstr filesep '..' filesep '..' filesep '..']; % Root of the Nemoh project
cwd = pwd();
search_locations = {nemoh_root; cwd};
nemoh_location = '';
for location=search_locations'
    match = dir([location{1} filesep '**' filesep 'mesh' executable_suffix]);
    for i=1:size(match,1)
        allMatch = [isfile([match(i).folder filesep 'mesh' executable_suffix])
                    isfile([match(i).folder filesep 'hydrosCal' executable_suffix])
                    isfile([match(i).folder filesep 'preProc' executable_suffix])
                    isfile([match(i).folder filesep 'solver' executable_suffix])
                    isfile([match(i).folder filesep 'postProc' executable_suffix])];
        if QTF
            allMatch = [allMatch
                        isfile([match(i).folder filesep 'QTFpreProc' executable_suffix])
                        isfile([match(i).folder filesep 'QTFsolver' executable_suffix])
                        isfile([match(i).folder filesep 'QTFpostProc' executable_suffix])];
        end
        if all(allMatch)
            nemoh_location = match(i).folder;
            disp(['Nemoh programs found in: ' nemoh_location]);
            break
        end
    end
    if ~isempty(nemoh_location)
        break
    end
end
if isfolder(nemoh_location)
    PATH = getenv('PATH');
    PATH = [PATH pathsep nemoh_location];
    setenv('PATH', PATH);
    disp(['Added directory ' nemoh_location ' to system PATH (for this session only).']);
    found = true;
    return
else
    disp('Nemoh programs not found in usual locations (working directory and Nemoh directories)')
end

% Prompt the user for a location and add it to the PATH
path = uigetdir('', "Nemoh programs directory");
if path ~= 0
    allMatch = [isfile([path filesep 'mesh' executable_suffix])
                isfile([path filesep 'hydrosCal' executable_suffix])
                isfile([path filesep 'preProc' executable_suffix])
                isfile([path filesep 'solver' executable_suffix])
                isfile([path filesep 'postProc' executable_suffix])];
    if QTF
        allMatch = [allMatch
                    isfile([path filesep 'QTFpreProc' executable_suffix])
                    isfile([path filesep 'QTFsolver' executable_suffix])
                    isfile([path filesep 'QTFpostProc' executable_suffix])];
    end
    if all(allMatch)
        nemoh_location = path;
        disp(['Nemoh programs found in: ' nemoh_location]);
        PATH = getenv('PATH');
        PATH = [PATH pathsep nemoh_location];
        setenv('PATH', PATH);
        disp(['Added directory ' nemoh_location ' to system PATH (for this session only).']);
        disp('Please consider permanently adding the Nemoh programs to your system PATH to avoid being prompted next time.')
        found = true;
        return
    else
        disp('The path you entered does not contain all required Nemoh programs.')
    end
else
    disp('Failed to find Nemoh programs.')
    found = false;
    return
end