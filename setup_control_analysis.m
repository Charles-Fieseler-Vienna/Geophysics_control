

%% Add current directory to path
what_am_i_setting_up = 'Control_analysis';
fprintf('Adding %s directory and subfolders to path\n',...
    what_am_i_setting_up)
folders_to_add = strsplit(genpath(pwd), pathsep);
% If we have git here, we'll have tons of folders we don't want
% Also don't want to add folders like 'scratch'
unwanted_folders = {'.git', 'scratch', 'pictures', 'doc', ...
    'intermediate_raw', 'figures'};
for i=1:length(unwanted_folders)
    unwanted_str = unwanted_folders{i};
    folders_to_add = folders_to_add(...
        cellfun(@(x)~contains(x,unwanted_str),folders_to_add));
end
% Clean up empty folder name
for i=1:length(folders_to_add)
    if isempty(folders_to_add{i})
        folders_to_add(i) = [];
    end
end

addpath(strjoin(folders_to_add, pathsep))
%==========================================================================


