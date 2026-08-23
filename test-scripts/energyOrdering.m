%% Signals to construct
% Case 1 - Closely spaced frequencies.
% NOTE: Cascade through increasing mode orders to observe how amplitude,
% decay factors, energies in the signal converge with appropriate model
% order selection.
clear;
clc;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;

f_true = [300, 300, 310, 500, 600, 1000, 2390];

amp_true = [1.5, 0.3, 0.75, 1, 1, 1, 1];

alpha_true = 50 : 10 : 110;
zeta_true = alpha_true ./ (f_true * 2 * pi);

phase_true = deg2rad(0 : 10 : 60);

M = length(f_true);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amp_true(m) * exp(-alpha_true(m) * t) .* ...
        sin(2 * pi * f_true(m) * t + phase_true(m));
    modeEnergyTrue(m) = sum(sigs(m, :) .^ 2);
end

modeIndex = sort(modeEnergyTrue, 'descend');
sig = sum(sigs);

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = true;
opts.filterOverdamped = false;
svdTolerance = 6;

eraLimit = 0.1;
fmin = 0;
fmax = 2500;

[~, minPhaseSig] = rceps(sig);

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, fEra, zetaEra, ampEra, phaseEra, modes_era, debug] = ...
    eigensystemRealisation(minPhaseSig, svdTolerance, Fs, opts);

poles = debug.Modal.poles;
residues = debug.Modal.residues;
alphaEra = debug.Modal.decayRates;

ii = find(fEra > fmin & fEra < fmax); % within analysis band
fEra = fEra(ii); ampEra = ampEra(ii); alphaEra = alphaEra(ii); phaseEra = phaseEra(ii);

ii = find(ampEra > eraLimit);
fEra = fEra(ii); alphaEra = alphaEra(ii); ampEra = ampEra(ii); phaseEra = phaseEra(ii);

phaseEra = phaseEra+pi/2;

inputImpulse = zeros(1, Ns);
inputImpulse(1) = 1;

% [B, A] = buildFilterCoefficients(poles, residues, sig(1));
[bEra, aEra] = generateIIMcoefficients(fEra, alphaEra, ampEra, phaseEra, dt);

[parallelFilterSignal, parallelEraModes] = parallelFilter(bEra, aEra, inputImpulse);
ERAmodeEnergy = sum(parallelEraModes .^ 2, 2);
[~, sortedEnergyIdx] = sort(ERAmodeEnergy, 'descend');
fEra = fEra(sortedEnergyIdx); alphaEra = alphaEra(sortedEnergyIdx);
ampEra = ampEra(sortedEnergyIdx); phaseEra = phaseEra(sortedEnergyIdx);

[freq, growth, amp, phase]=fdm_FAST(minPhaseSig(:)', Fs, fmin, fmax, 0.05);

ampE = 2*amp;      % correct for 2*cos(X) = exp(+j*X) + exp(-j*X)
alpE = -growth;       % flip growth to obtain attenuation
phaE = phase + 90;   % correct for sin rather than cos
ii = find(phaE > 180);  % ensure phase lying in [-180,180]
phaE(ii) = phaE(ii) - 360;

% REMOVING SPURIOUS COMPONENTS %%%%%
freES = freq; ampES = ampE; alpES = alpE; phaES = phaE;
ii = find(freES > fmin & freES < fmax); % within analysis band
freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii); 
ii = find(ampES > eraLimit);                % significant amplitude
freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii); 
ii = find(alpES > 0);                   % not growing
freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii);

nSigs = length(freES);
phaES = deg2rad(phaES);

[~, fdmSignal, ~, ~] = createImpulseResponse(nSigs, Ns, Fs, Frequency=freES, ...
    alpha=alpES, amplitude=ampES, phase=phaES, SNR=20);

[fdmFft, fdmFreqAxis] = singleSidedFft(fdmSignal, Fs);

% Generate parallel filter signals
[b, a] = generateIIMcoefficients(freES, alpES, ampES, phaES, dt);
[parallelFilterSignalFDM, parallelFDMModes] = parallelFilter(b, a, inputImpulse);
FDMmodeEnergy = sum(parallelFDMModes .^ 2, 2);
[~, sortedEnergyIdx] = sort(FDMmodeEnergy, 'descend');
freES = freES(sortedEnergyIdx); alpES = alpES(sortedEnergyIdx);
ampES = ampES(sortedEnergyIdx); phaES = phaES(sortedEnergyIdx);


% --- Match ERA modes to true modes by nearest frequency ---
% n = numel(f_true);
% idx = zeros(n,1);
% for k = 1:n
%     [~, idx(k)] = min(abs(debug.Reconstruction.modeEnergies - modeEnergyTrue(k)));
% end

modeIdx = 1 : 7;

f_e     = fEra(idx);
zeta_e  = zetaEra(idx);
amp_e   = ampEra(idx);
phase_e = phaseEra(idx) + pi/2;

% Errors
f_err     = round(100 * (f_e(:) - f_true(:)) ./ f_true(:), 2);
zeta_err  = round(100 * (zeta_e(:) - zeta_true(:)) ./ zeta_true(:), 2);
amp_err   = round(100 * (amp_e(:)  - amp_true(:))  ./ amp_true(:), 2);
phase_err = round(wrapToPi(phase_e(:) - phase_true(:)), 2);   % rad

T = table(modeIdx', ...
    f_true(:),     f_e(:),     f_err, ...
    zeta_true(:),  zeta_e(:),  zeta_err, ...
    amp_true(:),   amp_e(:),   amp_err, ...
    phase_true(:), phase_e(:), phase_err, ...
    'VariableNames', {'Mode', ...
    'Freq_true_Hz','Freq_ERA_Hz','Freq_err_pct', ...
    'Alpha_true','Alpha_ERA','Alpha_err_pct', ...
    'Amp_true','Amp_ERA','Amp_err_pct', ...
    'Phase_true_rad','Phase_ERA_rad','Phase_err_rad'});
 
format short g
newline
disp(T)

f = figure;
tl = tiledlayout(2, 2);
t1 = nexttile(tl);
plot(t, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original Signal')
hold on
plot(t, outputSignal, '--r', 'LineWidth', 1.3, 'DisplayName', 'ERA reconstructed signal')

t2 = nexttile(tl);
plot(t, sig, '--r', 'LineWidth', 1.5, 'DisplayName', 'ERA')
hold on
plot(t, parallelFilterSignalFDM, '-.b', 'LineWidth', 1.5, 'DisplayName', 'FDM')
xlabel(tl, 'Time (s)')
ylabel(tl, 'Amplitude')

f1 = figure;
tl = tiledlayout(1, 2);
nexttile
hold on

sv = debug.SVD.singularValues;
nSingularValues = length(sv);
svIndex = 1 : nSingularValues;
plot(svIndex, sv, '.k', 'DisplayName', 'Singular values')
plot(svIndex(1 : 14), sv(1 : 14), 'or', 'LineWidth',1.5)
ylabel('Singular value')
xlabel('Singular value index')
xlim([0 50])
nexttile
plot(svIndex, sv./sv(1), 'o', 'MarkerFaceColor','auto')
xlim([0 50])
yline(0.1, 'k--', "SVD Tolerance: " + string(0.1))
ylabel('Normalised singular value')
xlabel('Singular value index')

setFigureStyle()

%% Subfunctions
function [b, a] = generateIIMcoefficients(frequency, alpha, amplitude, ...
    phase, dt)

a0 = ones(length(frequency), 1);

R = exp(-alpha .* dt);
S = sin(2*pi*frequency.*dt);
C = cos(2*pi*frequency.*dt);

b0 = amplitude .* sin(phase);
b1 = amplitude .* R .* (cos(phase) .* S - sin(phase) .* C);
a1 = -2 .* R .* C;
a2 = R .^2;

b = [b0 b1];
a = [a0 a1 a2];

end

function setFigureStyle(varargin)
% SETFIGURESTYLE  Apply professional formatting to the current figure.
%
%   setFigureStyle() applies default professional settings.
%
%   setFigureStyle(Name, Value, ...) allows customisation via name-value pairs:
%
%   Parameters
%   ----------
%   'Width'       : Figure width in centimetres          (default: 16)
%   'Height'      : Figure height in centimetres         (default: 10)
%   'FontName'    : Font used for all text               (default: 'Arial')
%   'FontSize'    : Base font size in points             (default: 9)
%   'TitleSize'   : Title font size in points            (default: 10)
%   'LineWidth'   : Default line width for axes box      (default: 0.75)
%   'Style'       : Preset style string:
%                     'journal'  – single-column, minimal (default)
%                     'report'   – slightly larger, grid on
%                     'poster'   – large fonts, bold lines
%   'ColorOrder'  : N×3 RGB matrix of line colours.
%                   Defaults to a colourblind-friendly palette.
%   'Grid'        : Show grid lines – true | false       (default: false)
%   'Box'         : Show axes box    – true | false       (default: false)
%   'Target'      : Handle to a specific figure or axes  (default: gcf/gca)
%
%   Example
%   -------
%       x = linspace(0, 2*pi, 200);
%       figure; plot(x, sin(x), x, cos(x));
%       xlabel('x (rad)'); ylabel('Amplitude');
%       title('Sine and Cosine');
%       legend('sin(x)', 'cos(x)');
%       setFigureStyle('Style', 'journal', 'Width', 8.5, 'Height', 6);
 
% ── Parse inputs ──────────────────────────────────────────────────────────
p = inputParser;
p.CaseSensitive = false;
 
% Colourblind-friendly palette (Wong 2011, Nature Methods)
defaultColors = [
    0.00  0.45  0.70;   % blue
    0.84  0.37  0.00;   % vermillion
    0.00  0.62  0.45;   % green
    0.94  0.89  0.26;   % yellow
    0.34  0.71  0.91;   % sky blue
    0.90  0.62  0.00;   % orange
    0.80  0.47  0.65;   % reddish purple
];
 
addParameter(p, 'Width',      16,           @isnumeric);
addParameter(p, 'Height',     10,           @isnumeric);
addParameter(p, 'FontName',   'Arial',      @ischar);
addParameter(p, 'FontSize',   9,            @isnumeric);
addParameter(p, 'TitleSize',  10,           @isnumeric);
addParameter(p, 'LineWidth',  0.75,         @isnumeric);
addParameter(p, 'Style',      'journal',    @ischar);
addParameter(p, 'ColorOrder', defaultColors,@isnumeric);
addParameter(p, 'Grid',       false,        @islogical);
addParameter(p, 'Box',        false,        @islogical);
addParameter(p, 'Target',     [],           @(h) isempty(h) || ishghandle(h));
 
parse(p, varargin{:});
opts = p.Results;
 
% ── Apply style presets ───────────────────────────────────────────────────
switch lower(opts.Style)
    case 'journal'
        % Compact, single-column width, no grid
        opts.Width     = 8.5;
        opts.Height    = 6.0;
        opts.FontSize  = 8;
        opts.TitleSize = 9;
        opts.Grid      = false;
    case 'report'
        % Full-width with subtle grid
        opts.Width     = 16;
        opts.Height    = 10;
        opts.FontSize  = 9;
        opts.TitleSize = 10;
        opts.Grid      = true;
    case 'poster'
        % Large, bold for A0/A1 posters
        opts.Width     = 20;
        opts.Height    = 14;
        opts.FontSize  = 14;
        opts.TitleSize = 16;
        opts.LineWidth = 1.5;
        opts.Grid      = false;
    otherwise
        warning('setFigureStyle:unknownStyle', ...
            'Unknown style "%s". Using provided or default values.', opts.Style);
end
 
% ── Identify figure and axes ──────────────────────────────────────────────
if isempty(opts.Target)
    hFig = gcf;
    hAxes = findall(hFig, 'Type', 'axes');
elseif strcmp(get(opts.Target, 'Type'), 'figure')
    hFig  = opts.Target;
    hAxes = findall(hFig, 'Type', 'axes');
elseif strcmp(get(opts.Target, 'Type'), 'axes')
    hFig  = ancestor(opts.Target, 'figure');
    hAxes = opts.Target;
else
    error('setFigureStyle:invalidTarget', ...
        'Target must be a figure or axes handle.');
end
 
% ── Figure size (set in cm, convert to normalised paper units) ────────────
set(hFig, 'Units', 'centimeters');
pos = get(hFig, 'Position');
set(hFig, 'Position', [pos(1), pos(2), opts.Width, opts.Height]);
 
% Match paper size for PDF/EPS export
set(hFig, ...
    'PaperUnits',      'centimeters', ...
    'PaperSize',       [opts.Width, opts.Height], ...
    'PaperPosition',   [0, 0, opts.Width, opts.Height], ...
    'PaperPositionMode','manual');
 
% ── Apply settings to every axes in the figure ────────────────────────────
for i = 1:numel(hAxes)
    ax = hAxes(i);
 
    % Font
    set(ax, ...
        'FontName',  opts.FontName, ...
        'FontSize',  opts.FontSize);
 
    % Axes appearance
    set(ax, ...
        'LineWidth',   opts.LineWidth, ...
        'Box',         onoff(opts.Box), ...
        'XGrid',       onoff(opts.Grid), ...
        'YGrid',       onoff(opts.Grid), ...
        'ZGrid',       onoff(opts.Grid), ...
        'GridAlpha',   0.3, ...
        'TickDir',     'out', ...
        'TickLength',  [0.015, 0.025]);
 
    % Colour order for lines
    set(ax, 'ColorOrder', opts.ColorOrder, 'ColorOrderIndex', 1);
 
    % Labels and title font sizes
    set(get(ax, 'XLabel'), 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    set(get(ax, 'YLabel'), 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    set(get(ax, 'ZLabel'), 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    set(get(ax, 'Title'),  'FontName', opts.FontName, 'FontSize', opts.TitleSize, ...
        'FontWeight', 'bold');
 
    % Legend (if present)
    hLeg = legend(ax);
    if ~isempty(hLeg) && ishghandle(hLeg)
        set(hLeg, ...
            'FontName', opts.FontName, ...
            'FontSize', opts.FontSize, ...
            'Box',      'off', ...
            'Location', 'best');
    end
end
 
% ── Use renderer best suited for vector output ────────────────────────────
set(hFig, 'Renderer', 'painters');
 
end % setFigureStyle
 
% ── Helper: logical → 'on'/'off' ─────────────────────────────────────────
function s = onoff(tf)
    if tf; s = 'on'; else; s = 'off'; end
end

function exportFigure(filename, varargin)
% EXPORTFIGURE  Save the current figure to one or more file formats.
%
%   exportFigure(filename) saves the current figure as a 300 dpi PNG using
%   the filename stem provided (extension is added automatically).
%
%   exportFigure(filename, Name, Value, ...) allows customisation:
%
%   Parameters
%   ----------
%   'Formats'   : Cell array of format strings. Supported values:
%                   'png'  – raster, suitable for Word/PowerPoint reports
%                   'pdf'  – vector, suitable for LaTeX / print
%                   'eps'  – vector, high-quality journal submission
%                   'svg'  – vector, web / Inkscape editing
%                   'tif'  – raster TIFF, common for journal submissions
%                   'fig'  – MATLAB .fig source file
%                 Default: {'png'}
%   'DPI'       : Resolution in dots per inch for raster formats (default: 300)
%   'Target'    : Handle to a specific figure (default: gcf)
%   'OutDir'    : Output directory path (default: current directory)
%   'Overwrite' : Overwrite existing files – true | false (default: true)
%   'OpenDir'   : Open the output folder after export – true | false
%                 (default: false)
%
%   Returns
%   -------
%   A cell array of the full file paths written is printed to the command
%   window. The function also returns this list as an optional output.
%
%   Example
%   -------
%       figure; plot(rand(1,50));
%       setFigureStyle('Style', 'journal');
%       exportFigure('my_figure', ...
%           'Formats', {'png', 'pdf', 'eps'}, ...
%           'DPI',     600, ...
%           'OutDir',  'figures');
 
% ── Parse inputs ──────────────────────────────────────────────────────────
p = inputParser;
p.CaseSensitive = false;
 
addRequired( p, 'filename',  @ischar);
addParameter(p, 'Formats',   {'png'},  @iscell);
addParameter(p, 'DPI',       300,      @isnumeric);
addParameter(p, 'Target',    [],       @(h) isempty(h) || ishghandle(h));
addParameter(p, 'OutDir',    '.',      @ischar);
addParameter(p, 'Overwrite', true,     @islogical);
addParameter(p, 'OpenDir',   false,    @islogical);
 
parse(p, filename, varargin{:});
opts = p.Results;
 
% ── Resolve figure handle ─────────────────────────────────────────────────
if isempty(opts.Target)
    hFig = gcf;
else
    hFig = opts.Target;
end
 
if ~ishghandle(hFig) || ~strcmp(get(hFig, 'Type'), 'figure')
    error('exportFigure:invalidTarget', 'Target must be a valid figure handle.');
end
 
% ── Ensure output directory exists ────────────────────────────────────────
if ~isfolder(opts.OutDir)
    mkdir(opts.OutDir);
    fprintf('exportFigure: created output directory "%s"\n', opts.OutDir);
end
 
% Strip any extension the user may have included in filename
[~, stem, ~] = fileparts(opts.filename);
 
% ── Format definitions ────────────────────────────────────────────────────
% Maps format string → {print driver, file extension}
formatMap = struct( ...
    'png',  {{'-dpng',             'png'}}, ...
    'pdf',  {{'-dpdf',             'pdf'}}, ...
    'eps',  {{'-depsc',            'eps'}}, ...
    'svg',  {{'-dsvg',             'svg'}}, ...
    'tif',  {{'-dtiff',            'tif'}}, ...
    'fig',  {{'-dfig',             'fig'}}  ...  % handled separately
);
 
supportedFormats = fieldnames(formatMap);
 
% ── Export loop ───────────────────────────────────────────────────────────
writtenFiles = {};
 
for i = 1:numel(opts.Formats)
    fmt = lower(strtrim(opts.Formats{i}));
 
    % Validate format
    if ~ismember(fmt, supportedFormats)
        warning('exportFigure:unknownFormat', ...
            'Format "%s" is not supported and will be skipped.', fmt);
        continue
    end
 
    outPath = fullfile(opts.OutDir, [stem, '.', fmt]);
 
    % Overwrite check
    if ~opts.Overwrite && isfile(outPath)
        warning('exportFigure:fileExists', ...
            'File "%s" already exists and Overwrite is false. Skipping.', outPath);
        continue
    end
 
    fprintf('exportFigure: writing %s ... ', outPath);
 
    try
        if strcmp(fmt, 'fig')
            % savefig is the correct route for .fig files
            savefig(hFig, outPath);
 
        elseif ismember(fmt, {'png', 'tif'})
            % Raster: honour DPI setting
            dpiStr = sprintf('-r%d', opts.DPI);
            print(hFig, outPath, formatMap.(fmt){1}, dpiStr);
 
        else
            % Vector: DPI is ignored (use 0 to avoid artefacts on some versions)
            print(hFig, outPath, formatMap.(fmt){1}, '-r0');
        end
 
        fprintf('done.\n');
        writtenFiles{end+1} = outPath; %#ok<AGROW>
 
    catch ME
        fprintf('FAILED.\n');
        warning('exportFigure:writeFailed', ...
            'Could not write "%s": %s', outPath, ME.message);
    end
end
 
% ── Summary ───────────────────────────────────────────────────────────────
if isempty(writtenFiles)
    warning('exportFigure:noFilesWritten', 'No files were successfully exported.');
else
    fprintf('\nexportFigure: %d file(s) saved to "%s":\n', ...
        numel(writtenFiles), fullfile(pwd, opts.OutDir));
    for i = 1:numel(writtenFiles)
        fprintf('  %s\n', writtenFiles{i});
    end
end
 
% ── Optionally open output folder ─────────────────────────────────────────
if opts.OpenDir
    if ispc
        winopen(opts.OutDir);
    elseif ismac
        system(['open "', opts.OutDir, '"']);
    else
        system(['xdg-open "', opts.OutDir, '" &']);
    end
end
 
end % exportFigure
