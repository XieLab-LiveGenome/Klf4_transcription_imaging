%% concat_and_plot_cond_size.m
% Vertically concatenates the 4 condensate fields across multiple .mat files,
% saves a compiled .mat file, and plots the pooled diameter histogram.
% Also accumulates numNuclei across files for condensates-per-nucleus reporting.
clear; clc;
%% ---------------- User settings ----------------
folder  = '/Users/janaa/Desktop/MS2 transcription/Compiled data/Condensate size mat files /MED14/';
outName = 'concatenated_cond_size_filtered.mat';   % compiled output file
fields = {'condensate_feret','condensate_AR','condensate_FWHM','condensate_FWHM_r'};
%% ---------------- Gather files ----------------
files = dir(fullfile(folder, '*.mat'));
files = files(~strcmp({files.name}, outName));     % don't re-ingest the compiled file
assert(~isempty(files), 'No .mat files found in: %s', folder);
%% ---------------- Concatenate ----------------
for k = 1:numel(fields)
    C.(fields{k}) = [];
end
nFiles     = numel(files);
fileNames  = cell(nFiles,1);
fileMean   = nan(nFiles,1);
fileSD     = nan(nFiles,1);
fileMedian = nan(nFiles,1);
fileN      = zeros(nFiles,1);
fileNuclei = zeros(nFiles,1);
for i = 1:nFiles
    fpath = fullfile(files(i).folder, files(i).name);
    S = load(fpath);

    % pull numNuclei (sum in case it's stored per-FOV as a vector)
    if isfield(S, 'numNuclei')
        fileNuclei(i) = sum(S.numNuclei(:));
    else
        warning('Field "numNuclei" missing in %s', files(i).name);
        fileNuclei(i) = NaN;
    end

    for k = 1:numel(fields)
        f = fields{k};
        if isfield(S, f)
            C.(f) = [C.(f); S.(f)(:)];             % force column, then vertical concat
        else
            warning('Field "%s" missing in %s', f, files(i).name);
        end
    end
    % per-file stats on condensate diameter (Feret)
    d = S.condensate_feret(:);
    d = d(isfinite(d));
    fileNames{i}  = files(i).name;
    fileN(i)      = numel(d);
    fileMean(i)   = mean(d);
    fileSD(i)     = std(d);
    fileMedian(i) = median(d);
end
condensate_feret  = C.condensate_feret;
condensate_AR     = C.condensate_AR;
condensate_FWHM   = C.condensate_FWHM;
condensate_FWHM_r = C.condensate_FWHM_r;
totalNuclei       = sum(fileNuclei, 'omitnan');

save(fullfile(folder, outName), 'condensate_feret', 'condensate_AR', ...
    'condensate_FWHM', 'condensate_FWHM_r', 'totalNuclei', 'fileNuclei');
%% ---------------- Mean condensate size: per file + pooled ----------------
dAll        = condensate_feret(isfinite(condensate_feret));
pooledMean  = mean(dAll);
pooledSD    = std(dAll);
pooledMed   = median(dAll);
pooledN     = numel(dAll);
meanOfMeans = mean(fileMean);          % unweighted average across files
semOfMeans  = std(fileMean)/sqrt(nFiles);

fprintf('\n===== Mean condensate diameter (nm) =====\n');
fprintf('%-45s %8s %8s %10s %10s %10s %10s\n', ...
        'File', 'nNuc', 'n', 'Mean', 'SD', 'Median', 'Cond/Nuc');
fprintf('%s\n', repmat('-', 1, 105));
for i = 1:nFiles
    if fileNuclei(i) > 0
        cpn = fileN(i)/fileNuclei(i);
    else
        cpn = NaN;
    end
    fprintf('%-45s %8d %8d %10.1f %10.1f %10.1f %10.2f\n', ...
            fileNames{i}, fileNuclei(i), fileN(i), fileMean(i), fileSD(i), ...
            fileMedian(i), cpn);
end
fprintf('%s\n', repmat('-', 1, 105));
fprintf('%-45s %8d %8d %10.1f %10.1f %10.1f %10.2f\n', ...
        'POOLED (all condensates)', totalNuclei, pooledN, pooledMean, pooledSD, ...
        pooledMed, pooledN/totalNuclei);
fprintf('%-45s %8s %8d %10.1f %10.1f %10s %10s\n', ...
        'MEAN OF PER-FILE MEANS (+/- SEM)', '-', nFiles, meanOfMeans, semOfMeans, '-', '-');

fprintf('\nCompiled %d condensates from %d nuclei across %d files --> %s\n', ...
        pooledN, totalNuclei, nFiles, outName);

statsTable = table(fileNames, fileNuclei, fileN, fileMean, fileSD, fileMedian, ...
    'VariableNames', {'File','N_nuclei','N_condensates','Mean_nm','SD_nm','Median_nm'});

%% ---------------- Plot ----------------
ctrl_MED14    = condensate_feret;
ctrl_MED14_AR = condensate_AR;
figure
x_min = 0;
x_max = 1000;
s     = 25;
h_ctrl_MED14 = histogram(ctrl_MED14, 'LineWidth', 4.0, 'BinWidth', s, ...
    'FaceColor', [0.2 0.7 0.7]);
xlim([x_min x_max]);
xlabel('Condensate diameter (nm)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Counts', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('MED14 condensates (n = %d, %d nuclei)', numel(ctrl_MED14), totalNuclei), ...
    'FontSize', 20, 'FontWeight', 'bold');