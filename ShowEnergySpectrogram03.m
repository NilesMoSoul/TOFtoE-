%% Waterfall_Ridges_Publication.m
% Build a publication-grade "ridgeline/waterfall" plot from your Scan4.txt.
% - Row with "K.E ... ps" has delays, next row is units, data start after that.
% - First column = Energy [eV]; each subsequent column = signal at that delay.

%% ---------------- User controls ----------------
fname        = 'X:\rawdata\CFEL\CFEL-PAHs\TOF Spectrum\20240531\Scan 4 outupt\TOF0111\Scan4.txt';   % path to your file
OutName      = 'Scan4_ridges';          % base name for saved figures

NormalizePerTrace = false;              % true: normalize each column to max=1
BaselineShift     = 0;                  % vertical baseline for all ridges (z=BaselineShift)
YSpacing          = 5.0;                % fixed spacing between traces (arbitrary units)
EdgeWidth         = 0.6;                % outline width of each ridge
FaceAlphaFill     = 0.65;               % 0..1 transparency for fill
FigSizeInches     = [7.5 5.2];          % [W H] inches for export

% View / orientation
Azimuth     = -45;                      % user-adjustable viewing angle -35
Elevation   = 65;                       % user-adjustable viewing angle
XLim        = [];                       % e.g. [9 28]
ShowGrid    = 'off';                    % 'on' or 'off'
FontSize    = 11;

%% ---------------- Read & parse file (robust header find) ----------------
raw = fileread(fname);
raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
raw = strrep(raw, sprintf('\r'), sprintf('\n'));
lines = regexp(raw, '\n', 'split');

% drop empty/whitespace-only lines
lines = lines(~cellfun(@(s) isempty(s) || ~isempty(regexp(s,'^\s*$','once')), lines));

% find header line that contains delays (has "ps" or starts with K.E)
hdrIdx = [];
for i = 1:min(10, numel(lines))
    li = strtrim(lines{i});
    if contains(li, 'ps', 'IgnoreCase', true) || startsWith(li, 'K.E')
        hdrIdx = i; break;
    end
end
if isempty(hdrIdx)
    error('Could not find the header row containing delays (looked for "ps" / "K.E").');
end

% parse delays (strip first token "K.E", remove "ps", handle decimal commas)
hdr = strtrim(lines{hdrIdx});
hdr = regexprep(hdr, '\t', ' ');
hdr = regexprep(hdr, '\s+', ' ');
hdrNoPs = regexprep(hdr, 'ps', '', 'ignorecase');
hdrNoPs = regexprep(hdrNoPs, '^[^\s]+[\s]+', ''); % drop first token
hdrNoPs = strrep(hdrNoPs, ',', '.');
delayStrs = regexp(hdrNoPs, '\s+', 'split');
delays_ps = str2double(delayStrs);
delays_ps = delays_ps(~isnan(delays_ps));
nDel = numel(delays_ps);
if nDel == 0, error('No delays parsed from header row.'); end

% units row assumed directly after header
dataStartIdx = hdrIdx + 2;
if dataStartIdx > numel(lines)
    error('No data rows found after the units row.');
end

% numeric table
datLines = lines(dataStartIdx:end);
datLines = strrep(datLines, ',', '.'); % decimal commas -> dots

E  = nan(numel(datLines),1);
Z  = nan(numel(datLines), nDel);
rc = 0;
for i = 1:numel(datLines)
    s = strtrim(datLines{i}); if isempty(s), continue; end
    toks = regexp(s, '\s+', 'split');
    nums = str2double(toks); nums = nums(~isnan(nums));
    if isempty(nums), continue; end

    rc = rc + 1;
    E(rc) = nums(1);
    take = min(nDel, max(0, numel(nums)-1));
    zz = nan(1,nDel);
    if take > 0, zz(1:take) = nums(2:1+take); end
    Z(rc,:) = zz;
end
E = E(1:rc); Z = Z(1:rc,:);
bad = isnan(E) | all(isnan(Z),2);
E(bad) = []; Z(bad,:) = [];

% sort by energy
[Es, idx] = sort(E(:));
Zs = Z(idx,:);

% optional per-trace normalization
if NormalizePerTrace
    for j = 1:nDel
        m = max(abs(Zs(:,j)), [], 'omitnan');
        if isfinite(m) && m > 0, Zs(:,j) = Zs(:,j) / m; end
    end
end

%% ---------------- Visual styling helpers ----------------
% Pleasant, varied colors (desaturated a bit for “filled but faded” look)
Ncm = max(nDel, 7);
cm  = parula(Ncm);
% desaturate and soften
cm  = 0.8*cm + 0.2;                 % lift toward white
cIdx = round(linspace(1, Ncm, nDel));

% function to darken a color slightly for the edge
darken = @(c, f) max(0, min(1, c .* f));  % f<1 darkens

% y positions: fixed spacing (not the actual delay values)
yPos = (0:nDel-1) * YSpacing;

%% ---------------- Plot ridgeline waterfall ----------------
fig = figure('Color','w','Units','inches','Position',[1 1 FigSizeInches]);
ax  = axes('Parent', fig);
hold(ax,'on');

for j = 1:nDel
    x = Es(:);
    y = yPos(j) * ones(size(x));
    z = BaselineShift + Zs(:,j);

    % Create a closed patch (baseline back to start)
    xv = [x; flipud(x)];
    yv = [y; flipud(y)];
    zv = [zeros(size(z))+BaselineShift; flipud(z)];

    fc = cm(cIdx(j), :);
    ec = darken(fc, 0.65);

    p = patch(ax, xv, yv, zv, 1); %#ok<NASGU>
    set(p, 'FaceColor', fc, 'FaceAlpha', FaceAlphaFill, ...
              'EdgeColor', ec, 'LineWidth', EdgeWidth);

    % optional: thin crest line on top for crispness
    plot3(ax, x, y, z, 'Color', ec, 'LineWidth', EdgeWidth);
end

% Axes cosmetics
if ~isempty(XLim), xlim(ax, XLim); else, xlim(ax, [min(Es) max(Es)]); end
ylim(ax, [min(yPos)-0.5*YSpacing, max(yPos)+0.5*YSpacing]);

% Hide z-axis content
zlim(ax, 'auto');
ax.ZTick = [];
ax.ZColor = [1 1 1];     % effectively invisible
zlabel(ax, '');

% Clean frame: no box, optional grid
ax.Box = 'off';
grid(ax, ShowGrid);
ax.Layer = 'top';
ax.FontSize = FontSize;
ax.LineWidth = 0.9;

% Labels
xlabel(ax, 'K.E (eV)',          'FontWeight','bold');
ylabel(ax, 'Time Delay (ps)',   'FontWeight','bold');

% Put actual delay values as Y tick labels but at fixed positions
ax.YTick = yPos;
ax.YTickLabel = compose('%.3g ps', delays_ps);

% Nice view where x-axis faces the screen more—user can tweak Azimuth/Elevation above
view(ax, [Azimuth Elevation]);

% Light shading: subtle grid plane on floor only (optional)
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.12;

% Improve z baseline contrast a bit
set(ax, 'Clipping','on');

%% ---------------- Export ----------------
pngFile = sprintf('%s.png', OutName);
pdfFile = sprintf('%s.pdf', OutName);
exportgraphics(ax, pngFile, 'Resolution', 600, 'BackgroundColor','white');
exportgraphics(ax, pdfFile, 'ContentType','vector', 'BackgroundColor','white');
disp("Saved: " + string(pngFile));
disp("Saved: " + string(pdfFile));
