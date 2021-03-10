function [] = MakeFileStructure(path)
% MAKEFILESTRUCTURE: create folders and move images to folders as required
% for further analysis scripts
%
% Copyright (c) 2021 Dominik Schumacher, Silke Bergeler
%
% INPUT:
% - path: path to fluorescence and PH images

Fluor_exist = input('Data set contains fluorescence images? Y/N ','s');

if not(exist([path 'PH'],'dir'))
    mkdir(path, 'PH')
end
if Fluor_exist == 'Y' && not(exist([path 'Fluor'],'dir'))
    mkdir(path, 'Fluor')
end

name_PH = input('Enter characteristic expression for phase contrast images (e.g. PH or _ch00): ','s');

% list all tif files (extension can be captital letters or not) 
d = dir([path '*' name_PH '*.TIF']);

% if there are matches with an ending .TIF, move these files
if any(cellfun(@(x) ~isempty(regexp(x,['\w*' name_PH '\w*.TIF'],'once')),{d.name}))
    movefile([path '*' name_PH '*.TIF'],[path 'PH'])
end
% if there are matches with an ending .tif, move these files   
if any(cellfun(@(x) ~isempty(regexp(x,['\w*' name_PH '\w*.tif'],'once')),{d.name}))
    movefile([path '*' name_PH '*.tif'],[path 'PH'])
end

if Fluor_exist == 'Y'
    name_fluor = input('Enter characteristic expression for fluorescence images (e.g. Fluor or _ch01): ','s');
    
    d = dir([path '*' name_fluor '*.TIF']);
    
    % if there are matches with an ending .TIF, move these files
    if any(cellfun(@(x) ~isempty(regexp(x,['\w*' name_fluor '\w*.TIF'],'once')),{d.name}))
        movefile([path '*' name_fluor '*.TIF'],[path 'Fluor'])
    end
    % if there are matches with an ending .tif, move these files
    if any(cellfun(@(x) ~isempty(regexp(x,['\w*' name_fluor '\w*.tif'],'once')),{d.name}))
        movefile([path '*' name_fluor '*.tif'],[path 'Fluor'])
    end
end

%% save PH images to stack

PH_files = dir([path,'PH/*.tif']);

if contains(PH_files(1).name,"_RAW_ch0")
    PH_names = arrayfun(@(x) x.name, PH_files, 'UniformOutput',false);
    PH_paths = cellfun(@(x) [path,'PH/',x], PH_names, 'UniformOutput',false);
else
    for i = 1:length(PH_files)
        % split the name (name + format)
        name_splitted = strsplit(PH_files(i).name,'.');
        % select the name
        name_splitted = name_splitted{1};
        % test if the last character is a number
        name_digits = isstrprop(name_splitted,'digit');
        if ~name_digits(end)
            ph_dir_sep = PH_files(i); % save i-th entry separately
            PH_files(i) = []; % delete i-th entry
            continue; % there should be at most one entry without a number at the end
        end
    end
    PH_names = sort_nat(arrayfun(@(x) x.name, PH_files, 'UniformOutput',false)); % sort names of files in natural order
    
    if exist('ph_dir_sep','var')
        % add the name without a number as first entry
        C = {ph_dir_sep.name, PH_names};
        C = vertcat(C{:});
        PH_paths = cellfun(@(x) [path,'PH/',x], ...
            C, 'UniformOutput',false);
    else
        PH_paths = cellfun(@(x) [path,'PH/',x], PH_names, 'UniformOutput',false);
    end
end

image = [];
for i = 1:length(PH_paths)
    image = cat(3, image, imread(PH_paths{i})); % concatenate the images along the 3rd dimension to get a 3D image matrix
end

saveimagestack(image, [path 'PH-stack.tif'])

%% save TXR images to stack

if Fluor_exist == 'Y'
    % looks for .tif files in the Fluor folder. Contains only TXR pictures!
    TXR_files = dir([path,'Fluor/*.tif']);
    
    if contains(TXR_files(1).name,"_RAW_ch0")
        TXR_names = arrayfun(@(x) x.name, TXR_files, 'UniformOutput',false);
        TXR_paths = cellfun(@(x) [path,'Fluor/',x], TXR_names, 'UniformOutput',false);
    else
        for i = 1:length(TXR_files)
            % split the name (name + format)
            name_splitted = strsplit(TXR_files(i).name,'.');
            % select the name
            name_splitted = name_splitted{1};
            % test if the last character is a number
            name_digits = isstrprop(name_splitted,'digit');
            if ~name_digits(end)
                TXR_dir_sep = TXR_files(i); % save i-th entry separately
                TXR_files(i) = []; % delete i-th entry
                continue; % there should be at most one entry without a number at the end
            end
        end
        TXR_names = sort_nat(arrayfun(@(x) x.name, TXR_files, 'UniformOutput',false)); % sort names of files in natural order
        
        if exist('TXR_dir_sep','var')
            % add the name without a number as first entry
            C = {TXR_dir_sep.name, TXR_names};
            C = vertcat(C{:});
            TXR_paths = cellfun(@(x) [path,'Fluor/',x], ...
                C, 'UniformOutput',false);
        else
            TXR_paths = cellfun(@(x) [path,'Fluor/',x], TXR_names, 'UniformOutput',false);
        end
    end
    
    image = [];
    for i = 1:length(TXR_paths)
        image = cat(3, image, imread(TXR_paths{i})); % concatenate the images along the 3rd dimension to get a 3D image matrix
    end
    
    saveimagestack(image, [path 'TXR-stack.tif'])
end

end
