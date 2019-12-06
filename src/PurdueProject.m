classdef PurdueProject
% Various file paths and whatnot for the Purdue project
%
%
% INPUTS - none
%
% OUTPUTS - Object with preprocessing methods
%
% EXAMPLES
%
%   EXAMPLE1
%
%
%   EXAMPLE2
%
%
%
% Dependencies
%   .m files, .mat files, and MATLAB products required: 
%
%   See also: OTHER_FUNCTION_NAME
%
%
%
% Author: Charles Fieseler
% University of Washington, Dept. of Physics
% Email address: charles.fieseler@gmail.com
% Website: coming soon
% Created: 05-Dec-2019
%========================================
    
    properties (SetAccess = private)
        % Go straight to dropbox instead
%         dat_foldername = "../dat/ae_data/"
        dat_foldername = "C:/Users/charl/Dropbox/Kutz_SINDY/AE_Data/"
        dat_subfolders = {...
            "Distributed/Sample_490/", ...
            "Localized/Localized481_waveforms/", ...
            "Mortar/Mortar470_waveforms/", ...
            "Sandstone/Sandstone1/"}
        
        % For intermediate variables
        intermediate_foldername = "../intermediate_raw/"
        
        %
        presentation_foldername = "C:/Users/charl/Documents/Current_work/Presentations/Purdue/";
    end
    
    properties (Dependent)
        example_dat
        mortar_fnames
        localized_fnames
        distributed_fnames
    end
    
    methods % For dependent variables
        function out = get.example_dat(self)
            % Gets an example mortar data file
            fname = self.dat_foldername + self.dat_subfolders{3};
            tmp = dir(fname);
            out = readtable(fname + tmp(3).name);
        end
        
        function out = get.mortar_fnames(self)
            % Gets an example mortar data file
            fname = self.dat_foldername + self.dat_subfolders{3};
            out = self.get_fnames(fname);
        end
        
        function out = get.distributed_fnames(self)
            % Gets an example distributed data file
            fname = self.dat_foldername + self.dat_subfolders{1};
            out = self.get_fnames(fname);
        end
            
        function out = get.localized_fnames(self)
            % Gets an example localized data file
            fname = self.dat_foldername + self.dat_subfolders{2};
            out = self.get_fnames(fname);
        end
        
        function out = get_fnames(~, fname)
            tmp = dir(fname);
            tmp = tmp(3:end-2); % Last two are text files
            out = cell(size(tmp));
            for i = 1:length(tmp)
                out{i} = fname + tmp(i).name;
            end
            out = out(cellfun(@(x) contains(x, '.csv'), out));
        end
    end
    
    methods % For data pre-processing
        function [dat, kept_dat_ind]  =...
                filter_by_activity(~, fnames, min_activity_thresh, trim,...
                verbose)
            assert(iscell(fnames),...
                'Must pass cell array of filenames')
            if ~exist('min_activity_thresh', 'var') || ...
                    isempty(min_activity_thresh)
                min_activity_thresh = 0.015;
            end
            if ~exist('trim', 'var')
                trim = 500;
            end
            if ~exist('verbose', 'var')
                verbose = true;
            end
            if verbose
                disp('Filtering by activity...')
            end
            num_datasets = length(fnames);
            kept_dat_ind = false(num_datasets, 1);
            dat = cell(num_datasets, 1);
            all_times = dat;
            for i = 1:num_datasets
                if verbose
                    fprintf('File %d/%d\n', i, num_datasets)
                end
                dat{i} = readtable(fnames{i});
                all_times{i} = dat{i}{1:end-trim, 1}';
                dat{i} = dat{i}{1:end-trim,2}';
                kept_dat_ind(i) = (max(dat{i}) > min_activity_thresh);
            end
            if verbose
                disp('Finished filtering')
            end

            dat = dat(kept_dat_ind);
        end
        
        function [dat, kept_dat_ind] = ...
                filter_active_sensing(~, fnames, max_activity_thresh)
            % Events that have activity near the maximum level of 10 are
            % actually "Active Sensing Tests" (ASTs), as opposed to
            % spontaneous events
            assert(iscell(fnames),...
                'Must pass cell array of filenames')
            if ~exist('min_activity_thresh', 'var') || ...
                    isempty(max_activity_thresh)
                max_activity_thresh = 9.9;
            end
            num_datasets = length(fnames);
            kept_dat_ind = false(num_datasets, 1);
            dat = cell(num_datasets, 1);
            for i = 1:num_datasets
                dat{i} = readtable(fnames{i});
                dat{i} = dat{i}{1:end,2}';
                kept_dat_ind(i) = (max(dat{i}) < max_activity_thresh);
            end

            dat = dat(kept_dat_ind);
        end
        
    end
    
    methods % For processing metadata
        function [out] = import_metadata_file(~, which_dataset)
            if strcmp(which_dataset, 'mortar')
                fname = 'C:\Users\charl\Documents\MATLAB\Collaborations\Purdue_analysis_git\dat\ae_data\mortar\Mortar470_waveforms\180225152101_events';
            else
                error('Unrecognized metadata')
            end
            out = readtable(fname);
            
            time_col = days(out{:,2}) + out{:,3};
            out = [out(:, 4:end) table(time_col)];
        end
        
        function ind = fname2metadataInd(~, fname, mdata)
            if contains(fname, '/')
                tmp = strsplit(fname, '/');
                fname = tmp{end};
            end
            if contains(fname, '.csv')
                fname = fname(1:end-4);
            end
            vals = strsplit(fname, '_');
            
            channel = str2num(vals{3});
            ind_in_channel = str2num(vals{4});
            
            all_ind = find(mdata{:,'CH'}==channel);
            ind = all_ind(ind_in_channel);
        end
        
        function m = fname2metadata(self, fname, mdata)
            ind = self.fname2metadataInd(fname, mdata);
            m = mdata(ind,:);
        end
        
        function t = read_time_of_test(~, fname)
            fid = fopen(fname{1});
            C = textscan(fid, '%s%f', 1, ...
                'Headerlines', 10, 'Delimiter', ':');
            fclose(fid);
            t = C{2};
        end
    end
end

