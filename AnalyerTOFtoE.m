%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOF_Streaking_Option1.m
% Generates TOF matrix and energy-calibrated streaking spectrum
% from a folder of .dat files (each with [TOF(s), Amplitude]).
%
% Only "Option 1" logic kept. No menus or slicing branch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear; tic

%% ---------------- Calibration block (select ONE set) --------------
% Example: Acceleration 1 V @ 4 kS/s
Uret = -0.141;
t0   = 13.327;
cCal = 733.193;

% --- Keep other calibrations here as commented presets for reference ---
% % Accel 0 V @ 10 kS/s (Std lenses)
% Uret = 0.68; t0 = 20.93; cCal = 670.75;

% % Accel 15.3 V
% Uret = 10.03; t0 = 27.66; cCal = 629.7;

% % Accel 5 V (OSLS)
% Uret = 3.482125; t0 = 22.309082; cCal = 679.848597;

%% --------------- Optional external data (unused but kept) ----------
% Safely load if present (won't error if missing)
if exist('tmap.dat','file'),      tmap = load('tmap.dat');       else, tmap = []; end
if exist('output_data.dat','file'), harmonics = load('output_data.dat'); else, harmonics = []; end

%% ----------------------- Directories -------------------------------
uiwait(msgbox('Please select the directory of the raw data (.dat files).','Raw Data Directory','help'));
rawDataDir = uigetdir();
if isequal(rawDataDir,0), error('Directory selection canceled.'); end
disp(['Raw data directory: ', rawDataDir]);

uiwait(msgbox('Please select the directory for storing the output files.','Output Directory','help'));
outputDir = uigetdir();
if isequal(outputDir,0), error('Directory selection canceled.'); end
disp(['Output directory: ', outputDir]);

%% ----------------------- Gather files ------------------------------
fileList = dir(fullfile(rawDataDir,'*.dat'));
if isempty(fileList), error('No .dat files found in: %s', rawDataDir); end

% Natural sort by filename
fileNames = {fileList.name};
fileNames = natsort_cell(fileNames);
fileList  = fileList(ismember({fileList.name}, fileNames));
[~, order] = ismember(fileNames, {fileList.name}); %#ok<ASGLU> % kept for clarity

numFiles = numel(fileList);
fprintf('Number of points: %d\n', numFiles);

%% --------------------- User inputs ---------------------------------
% Temporal scan step (fs) between successive spectra
Step_Size_fs = input('Please enter the step size between spectra (in femtoseconds): ');
if ~isscalar(Step_Size_fs) || ~isfinite(Step_Size_fs) || Step_Size_fs<=0
    error('Step size must be a positive number (fs).');
end

% TOF window (ns) to use for energy conversion & plotting
prompt = {'TOF window start (ns):','TOF window end (ns):'};
answ   = inputdlg(prompt,'TOF window for analysis',1,{'100','360'});
if isempty(answ), error('Canceled.'); end
TOF_ns_min = str2double(answ{1});
TOF_ns_max = str2double(answ{2});
if any(~isfinite([TOF_ns_min, TOF_ns_max])) || TOF_ns_min>=TOF_ns_max
    error('Invalid TOF window.');
end

%% ---------------- Build TOF matrix (amplitude vs file index) -------
dataCols = cell(numFiles,1);

for i = 1:numFiles
    fpath = fullfile(rawDataDir, fileNames{i});
    M = load(fpath);
    if size(M,2) < 2
        error('File %s does not have at least 2 columns.', fileNames{i});
    end
    amp = -M(:,2);                % negate as in original code
    amp = amp - min(amp);         % shift to positive
    dataCols{i} = amp;
end

maxRows = max(cellfun(@numel, dataCols));
TOFMatrix = NaN(maxRows, numFiles);
for i = 1:numFiles
    n = numel(dataCols{i});
    TOFMatrix(1:n,i) = dataCols{i};
end

% Save TOF matrix
tofMatPath = fullfile(outputDir,'TOFmatrix.dat');
writematrix(TOFMatrix, tofMatPath, 'Delimiter','\t');
disp(['TOF Matrix saved: ', tofMatPath]);

%% ---------------- Time axes ----------------------------------------
NP             = numFiles;                         % number of scanned points
scanRange_ps   = (Step_Size_fs * NP)/1000;         % total scan range in ps
scanStep_ps    = scanRange_ps / NP;                % step in ps
tau_ps         = 0:scanStep_ps:(scanRange_ps - scanStep_ps);  % length NP

% Read first file to get TOF sampling axis
firstPath      = fullfile(rawDataDir, fileNames{1});
firstData      = load(firstPath);
tof_ns_full    = firstData(:,1)*1e9;               % convert s -> ns
if ~issorted(tof_ns_full), [tof_ns_full, idxsrt] = sort(tof_ns_full); else, idxsrt = []; end

% Define window indices
start_row = find(tof_ns_full >= TOF_ns_min, 1, 'first');
end_row   = find(tof_ns_full <= TOF_ns_max, 1, 'last');
if isempty(start_row) || isempty(end_row) || end_row<=start_row
    error('Chosen TOF window not found in sampling axis.');
end

% Copy the first scanned spectrum to output for reference
copyfile(firstPath, fullfile(outputDir,'First Scanned Spectrum.dat'));

%% ---------------- Plot raw TOF “streaking trace” -------------------
figure('Color','w');
imagesc(tau_ps, tof_ns_full, TOFMatrix);
set(gca,'YDir','normal'); colormap jet; colorbar
title('TOF Photoelectron Spectrum');
xlabel('Time Delay (ps)'); ylabel('TOF (ns)');
ylim([TOF_ns_min TOF_ns_max]);

% Unique filename in outputDir
base1 = 'TOF Streaking Trace';
k=0; fn1 = fullfile(outputDir, sprintf('%s %03d.png',base1,k));
while exist(fn1,'file'), k=k+1; fn1 = fullfile(outputDir, sprintf('%s %03d.png',base1,k)); end
saveas(gcf, fn1);
fprintf('Plot saved: %s\n', fn1);

%% -------------- Energy calibration & spectrum matrix ---------------
% Build energy axis (from selected TOF window) using calibration:
%   E(eV) = (cCal^2) / ( (TOF_ns - t0)^2 ) - Uret
% Use amplitude Jacobian factor from your code:
%   H = 0.5 * A * cCal / ( (E+Uret)^(3/2) )

W_eV = (cCal^2) ./ ((tof_ns_full(start_row:end_row) - t0).^2) - Uret;  % energy axis (ascending with decreasing TOF)
W_eV = flip(W_eV);                                                     % flip to increasing energy upwards

EnergyCols = cell(numFiles,1);

for i = 1:numFiles
    fpath = fullfile(rawDataDir, fileNames{i});
    M = load(fpath);
    TOF_ns  = M(:,1)*1e9;
    Amp_raw = -M(:,2);

    % If first file was sorted, apply same sorting to be consistent
    if ~isempty(idxsrt)
        TOF_ns  = TOF_ns(idxsrt);
        Amp_raw = Amp_raw(idxsrt);
    end

    E_eV = (cCal^2) ./ ((TOF_ns(start_row:end_row) - t0).^2) - Uret;
    Hamp = 0.5 * Amp_raw(start_row:end_row) .* (cCal ./ ((E_eV + Uret).^(3/2)));

    % Clean and flip to match W_eV ordering
    valid = isfinite(E_eV) & isfinite(Hamp) & (E_eV > 0);
    E_eV  = E_eV(valid);
    Hamp  = Hamp(valid);

    EnergyCols{i} = flip(Hamp); % flip so rows correspond to ascending W_eV
end

maxRowsE = max(cellfun(@numel, EnergyCols));
EnergyMatrix = NaN(maxRowsE, numFiles);
for i = 1:numFiles
    n = numel(EnergyCols{i});
    EnergyMatrix(1:n, i) = EnergyCols{i};
end

% Save energy-calibrated matrix
engMatPath = fullfile(outputDir,'EnergyCalibratedMatrix.dat');
writematrix(EnergyMatrix, engMatPath, 'Delimiter','\t');
disp(['Energy Calibration Matrix saved: ', engMatPath]);

%% -------------- Plot photoelectron streaking spectrum --------------
figure('Color','w');
pcolor(tau_ps, W_eV, EnergyMatrix); shading interp
colormap jet; colorbar
title('Photoelectron Streaking Spectrum');
xlabel('Time Delay (ps)'); ylabel('K.E. (eV)');
ylim([5.2 25]);

base2 = 'Photoelectron Streaking Spectrum';
k=0; fn2 = fullfile(outputDir, sprintf('%s %03d.png',base2,k));
while exist(fn2,'file'), k=k+1; fn2 = fullfile(outputDir, sprintf('%s %03d.png',base2,k)); end
saveas(gcf, fn2);
fprintf('Plot saved: %s\n', fn2);

%% ---------------- Save scan parameters -----------------------------
paramNames  = {'Number of Scanned Points','Scan Range (ps)','Step Size (fs)','Time Step (ps)','TOF first row idx','TOF last row idx','TOF min (ns)','TOF max (ns)'};
paramValues = [NP, scanRange_ps, Step_Size_fs, scanStep_ps, start_row, end_row, TOF_ns_min, TOF_ns_max];

% Human-readable (name + value)
scanParamsPath = fullfile(outputDir,'ScanParameters.dat');
fid = fopen(scanParamsPath,'w');
for i=1:numel(paramNames)
    fprintf(fid,'%s\t%g\n', paramNames{i}, paramValues(i));
end
fclose(fid);
disp(['Parameters saved: ', scanParamsPath]);

% Machine-friendly (values only, one per line)
scanParamsValsPath = fullfile(outputDir,'ScanParameters1.dat');
fid = fopen(scanParamsValsPath,'w');
fprintf(fid,'%g\n', paramValues);
fclose(fid);
disp(['Parameters saved: ', scanParamsValsPath]);

%% ---------------- Save workspace snapshot --------------------------
save(fullfile(outputDir,'TOF_Streaking_Workspace.mat'), ...
    'tau_ps','tof_ns_full','W_eV','TOFMatrix','EnergyMatrix', ...
    'NP','Step_Size_fs','scanStep_ps','start_row','end_row', ...
    'Uret','t0','cCal','tmap','harmonics');

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- Helper: natural sort of cell array of filenames ----------
function C = natsort_cell(C)
% Simple natural sort: splits digit runs and compares numerically
tokens = cellfun(@(s) regexp(s,'(\d+)|(\D+)','match'), C, 'UniformOutput', false);
numtokens = cellfun(@numel, tokens);
maxT = max(numtokens);
key = zeros(numel(C), maxT*2); % [isNum val] pairs

for i = 1:numel(C)
    tks = tokens{i};
    row = [];
    for k = 1:numel(tks)
        tk = tks{k};
        isnum = ~isempty(regexp(tk,'^\d+$','once'));
        if isnum
            row = [row, 1, str2double(tk)]; %#ok<AGROW>
        else
            % map text token to rank via lowercased char codes packed into base-256
            lc = double(lower(tk));
            rank = 0;
            for ii=1:numel(lc)
                rank = rank*256 + lc(ii);
            end
            row = [row, 0, rank]; %#ok<AGROW>
        end
    end
    key(i,1:numel(row)) = row;
end

[~,idx] = sortrows(key);
C = C(idx);
end
