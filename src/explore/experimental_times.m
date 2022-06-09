%% Get the raw data (experimental times for Chven)

%% Get the indices of the fnames within the list
pp = PurdueProject();

out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig9_indices.mat";
ind_struct = load(out_fname);

%% Load the actual data: 0

fname = "mortar_fnames_all_raw_data";
dat0 = load(pp.intermediate_foldername + fname + ".mat");

for i = ind_struct.example_ind_original
    fname = dat0.fnames{i};
    open(fname)
    pause
end

%% Load the actual data: 1

fname = "mortar_new1_fnames_all_raw_data";
dat0 = load(pp.intermediate_foldername + fname + ".mat");

for i = ind_struct.example_ind_mortar1
    fname = dat0.fnames{i};
    open(fname)
end

%% Load the actual data: 2

fname = "mortar_new2_fnames_all_raw_data";
dat0 = load(pp.intermediate_foldername + fname + ".mat");

for i = ind_struct.example_ind_mortar2
    fname = dat0.fnames{i};
    open(fname)
end

%% Load the actual data: 3

fname = "mortar_new3_fnames_all_raw_data";
dat0 = load(pp.intermediate_foldername + fname + ".mat");

for i = ind_struct.example_ind_mortar3
    fname = dat0.fnames{i};
    open(fname)
end

% Plot one specifically
% dat = open(dat0.fnames{ind_struct.example_ind_mortar3(1)});
% plot(dat)