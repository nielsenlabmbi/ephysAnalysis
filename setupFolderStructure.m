function setupFolderStructure(rootFolder,createRaw,createOut)
% creates the default folder structure in the rootFolder.
% if rootFolder is not supplied, it creates the structure inside the
% current directory.
% err and settings folders are always created (in the current directory).
% createRaw and createOut default to 1 if not supplied.
% createRaw = 1 for local, createRaw = 2 for synology
% it is assumed that the OutShare folder structure in the synology will be
% manually created.

if ispc
    slash = '\';
else
    slash = '/';
end

if ~exist('createOut','var');           createOut = 1;      end     
if ~exist('createRaw','var');           createRaw = 1;      end    
if ~exist('rootFolder','var')
    rootFolder = '';    
elseif ~isempty(rootFolder) && ~strcmp(rootFolder(end),slash)
    rootFolder(end+1) = slash;
end    

if ~exist('settings','dir');    mkdir('settings');  end;
if ~exist('err','dir');         mkdir('err');       end;

if createOut
    if ~exist([rootFolder 'out'],'dir');                       mkdir([rootFolder 'out']);                     end;
    if ~exist([rootFolder 'out' slash 'dataStage'],'dir');     mkdir([rootFolder 'out' slash 'dataStage']);   end;
    if ~exist([rootFolder 'out' slash 'spikeStage'],'dir');    mkdir([rootFolder 'out' slash 'spikeStage']);  end;
    if ~exist([rootFolder 'out' slash 'summary'],'dir');       mkdir([rootFolder 'out' slash 'summary']);     end;
    if ~exist([rootFolder 'out' slash 'plots'],'dir');         mkdir([rootFolder 'out' slash 'plots']);       end;
end

if createRaw
    if ~exist([rootFolder 'raw'],'dir');                       mkdir([rootFolder 'raw']);                     end;
    if ~exist([rootFolder 'raw' slash 'data'],'dir');          mkdir([rootFolder 'raw' slash 'data']);        end;
    if ~exist([rootFolder 'raw' slash 'analyzer'],'dir');      mkdir([rootFolder 'raw' slash 'analyzer']);    end;
end