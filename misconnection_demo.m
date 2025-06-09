%% Install drEEM Toolbox
clc;
clear all;

% installation
% change this path to your drEEM tool box path
cd('/Users/caroline/Desktop/EEMs/drEEM-0.6.2');
dreeminstall;

%% Clear Workspace and Load Dataset
clearvars;            % Clear all variables
close all force;      % Close all figures

%% Preprocess EEM Data and Remove Scattering

% Navigate to the dataset directory
% change this path to your own dataset path
cd('/Users/caroline/Desktop/EEMs_DATASET/format_DATA/A');

% Load excitation and emission wavelength vectors
Ex = xlsread('A-01.xlsx', 'B37:FU37');     % Excitation wavelengths
Em = xlsread('A-01.xlsx', 'A38:A213');     % Emission wavelengths

% Read EEM dataset (multiple samples)
% Syntax: readineems(type, format, range, headers, display, outdata)
[X, Emmat, Exmat, flist, outdata] = readineems(1, 'A-*.xlsx', 'A37:FU213', [1 1], 0, 0);

% Assemble data into a drEEM-compatible 1x1 struct
% Syntax: assembledataset(X, Ex, Em, units, name1, data1, ...)
DS0 = assembledataset(X, Ex(1,:), Em(:,1), 'AU', 'flist', flist, []);

% Quick visual review of the dataset
eemreview(DS0);

% Crop data by trimming wavelength ranges
% Syntax: subdataset(dataset, outSample, outEm, outEx, backups)
DS1 = subdataset(DS0, [], DS0.Em > 505, DS0.Ex > 400);    % First trim
DS  = subdataset(DS1, [], DS0.Em > 500, DS0.Ex < 204);    % Second trim

% Review the trimmed dataset
eemreview(DS);

% Remove scattering components (may require multiple adjustments)
% Syntax: smootheem(data, Ray1, Ram1, Ray2, Ram2, NaNfilter, d2zero, freq, plot)
Xs = smootheem(DS, [12 5], [20 10], [15 15], [7 8], [1 1 1 1], [], 3382, 0);

% Visual check of scatter removal
eemreview(Xs);               % Close-up view of corrected regions
eemview(Xs, 'X', [8 5]);     % Overview plot if needed
