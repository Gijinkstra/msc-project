%% Manage interpreters
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
%% Signals to construct
% Case 1 - Closely spaced frequencies.
% NOTE: Cascade through increasing mode orders to observe how amplitude,
% decay factors, energies in the signal converge with appropriate model
% order selection.
clear;
clc;

fs = 10000;
dt = 1/fs;
nSamples = 1000;

f_true = [300, 300, 310, 500, 600, 1000, 2390];

amp_true = [1.5, 0.3, 0.75, 1, 1, 1, 1];

alpha_true = 50 : 10 : 110;
zeta_true = alpha_true ./ (f_true * 2 * pi);

phase_true = deg2rad(0 : 20 : 120);

nSignals = numel(f_true);

[allSignals, outputSignal, timeVector, debug] = ...
    createImpulseResponse(nSignals, nSamples, fs, SNR=0, Frequency=f_true, ...
    Alpha=alpha_true, Amplitude=amp_true, Phase=phase_true);

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = false;
svdTolerance = 7;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[eraSignal, fEra, dfEra, ampEra, phaseEra, modesEra, debug] = ...
    eigensystemRealisation(outputSignal, svdTolerance, fs, opts);

poles = debug.Modal.poles;
residues = debug.Modal.residues;

inputImpulse = zeros(1, nSamples);
inputImpulse(1) = 1;

[B, A] = buildFilterCoefficients(poles, residues, outputSignal(1));
parallelFilterSignal = parallelFilter(B, A, inputImpulse);

f1 = figure;
tl = tiledlayout(2, 1);
xlabel(tl, 'Time (s)', 'Interpreter', 'Latex')
ylabel(tl, 'Amplitude', 'Interpreter', 'Latex')

t1 = nexttile;
hold(t1, "on")

t2 = nexttile;
hold(t2, "on")

parallelFilterNMSE = sum((outputSignal(:) - ...
    parallelFilterSignal(:)) .^ 2) / sum(outputSignal(:) .^2);

% t1 — True vs ERA
plot(t1, timeVector, eraSignal,    '-',  'Color', [0.84 0.37 0.00 0.55], ...
     'LineWidth', 3, 'DisplayName', 'ERA signal')
plot(t1, timeVector, outputSignal, '-',  'Color', [0.00 0.45 0.70], ...
     'LineWidth', 1.2, 'DisplayName', 'True signal')
legend(t1, 'Location', 'southeast', 'Box', 'off')
% text(0.8, 0.9, "Reconstruction error: " + string(debug.Reconstruction.impulseNMSE))

% t2 — ERA vs Parallel Filter
plot(t2, timeVector, parallelFilterSignal, '-',  'Color', [0.00 0.62 0.45 0.55], ...
     'LineWidth', 3, 'DisplayName', 'Parallel Filter signal')
plot(t2, timeVector, eraSignal,            '-',  'Color', [0.84 0.37 0.00], ...
     'LineWidth', 1.2, 'DisplayName', 'ERA signal')
legend(t2, 'Location', 'southeast', 'Box', 'off')

modeIdx = 1:7;

phaseEra = phaseEra+pi/2;

% phase_true = round(phase_true, 2, "significant");
% alpha_true = round(alpha_true, 2, "significant");
% amp_true = round(amp_true, 2, "significant");
% f_true = round(f_true, 2, "significant");
% 
% phaseEra = round(phaseEra, 2, "significant");
% alpha_true = round(alpha_true, 2, "significant");
% amp_true = round(amp_true, 2, "significant");
% f_true = round(f_true, 2, "significant");

phaseEra([2, 1]) = phaseEra([1, 2]);
fEra([2, 1]) = fEra([1, 2]);
ampEra([2, 1]) = ampEra([1, 2]);
dfEra([2, 1]) = dfEra([1, 2]);

alpEra = dfEra .* (2*pi.*fEra);

% --- percentage errors (element-wise, signed) ---
% guard against divide-by-zero where a true value is 0 (e.g. phase mode 1)
pctErr = @(est, tru) 100 * (est - tru) ./ tru;

err_f     = pctErr(fEra(:),      f_true(:));
err_alpha = pctErr(alpEra(:),    alpha_true(:));
err_amp   = pctErr(ampEra(:),    amp_true(:));
err_phase = pctErr(phaseEra(:),  phase_true(:));

% --- assemble and print ---
Terr = table(modeIdx(:), f_true(:), err_f, err_alpha, err_amp, err_phase, ...
    'VariableNames', {'Mode','f_Hz','err_f','err_alpha','err_amp','err_phase'});
disp(Terr)

% --- emit LaTeX body rows (tab/ampersand separated, 3 sig figs on errors) ---
fprintf('\n--- LaTeX rows ---\n');
for k = 1:numel(modeIdx)
    fprintf('%d & %g & %.3g & %.3g & %.3g & %.3g \\\\\n', ...
        modeIdx(k), f_true(k), err_f(k), err_alpha(k), err_amp(k), err_phase(k));
end

% T = table(modeIdx', ...
%     f_true(:),     fEra(:), ...
%     alpha_true(:),  alpEra(:), ...
%     amp_true(:),   ampEra(:), ...
%     phase_true(:), phaseEra(:), ...
%     'VariableNames', {'Mode', ...
%     'True Frequency','ERA Frequency', ...
%     'True Alpha','ERA Alpha', ...
%     'True Amp','Era Amp', ...
%     'True Phase','ERA Phase'});
% 
% format short g
% newline
% disp(T)

setFigureStyle('style', 'Report')

f2 = figure;
til2 = tiledlayout(f2, 1, 2);

t1 = nexttile(til2);
hold(t1, "on")
xlabel("Frequency (Hz)")
ylabel("Decay")

t2 = nexttile(til2);
hold(t2, "on")
xlabel('Phase (rad)')
ylabel('Amplitude')

plot(t1, f_true, alpha_true, 'ok', 'LineWidth', 3, 'MarkerSize', 10, 'DisplayName', 'Simulated values')
plot(t1, fEra, alpEra, 'rx', 'LineWidth', 3,  'MarkerSize', 8, 'DisplayName', 'Simulated values')
legend

plot(t2, phase_true, amp_true, 'ok', 'LineWidth', 3, 'MarkerSize', 10, 'DisplayName', 'Simulated values')
plot(t2, phaseEra, ampEra, 'rx', 'LineWidth', 3, 'MarkerSize', 8, 'DisplayName', 'Simulated values')
legend

%% Noisy signal simulation
clear;
clc;

fs = 10000;
dt = 1/fs;
nSamples = 1000;
nSignals = 20;
SNR = 30;

rng(42, "twister");

[allSignals, outputSignal, timeVector, debugImpulse] = ...
    createImpulseResponse(nSignals, nSamples, fs, SNR=SNR, Debug=true);
[~, sortIdx] = sort(debugImpulse.frequency);

f_true = debugImpulse.frequency(sortIdx);
amp_true = debugImpulse.amplitude(sortIdx);
alpha_true = debugImpulse.alpha(sortIdx);
phase_true = debugImpulse.phase(sortIdx);

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = false;
svdTolerance = 20;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[eraSignal, fEra, dfEra, ampEra, phaseEra, modesEra, debug] = ...
    eigensystemRealisation(outputSignal, svdTolerance, fs, opts);

poles = debug.Modal.poles;
residues = debug.Modal.residues;

inputImpulse = zeros(1, nSamples);
inputImpulse(1) = 1;

[B, A] = buildFilterCoefficients(poles, residues, outputSignal(1));
parallelFilterSignal = parallelFilter(B, A, inputImpulse);

f1 = figure;
tl = tiledlayout(2, 1);
xlabel(tl, 'Time (s)', 'Interpreter', 'Latex')
ylabel(tl, 'Amplitude', 'Interpreter', 'Latex')

t1 = nexttile;
hold(t1, "on")

t2 = nexttile;
hold(t2, "on")

parallelFilterNMSE = sum((outputSignal(:) - ...
    parallelFilterSignal(:)) .^ 2) / sum(outputSignal(:) .^2);

% t1 — True vs ERA
plot(t1, timeVector, eraSignal,    '-',  'Color', [0.84 0.37 0.00 0.55], ...
     'LineWidth', 3, 'DisplayName', 'ERA signal')
plot(t1, timeVector, outputSignal, '-',  'Color', [0.00 0.45 0.70], ...
     'LineWidth', 1.2, 'DisplayName', 'True signal')
legend(t1, 'Location', 'southeast', 'Box', 'off')
% text(0.8, 0.9, "Reconstruction error: " + string(debug.Reconstruction.impulseNMSE))

% t2 — ERA vs Parallel Filter
plot(t2, timeVector, parallelFilterSignal, '-',  'Color', [0.00 0.62 0.45 0.55], ...
     'LineWidth', 3, 'DisplayName', 'Parallel Filter signal')
plot(t2, timeVector, eraSignal,            '-',  'Color', [0.84 0.37 0.00], ...
     'LineWidth', 1.2, 'DisplayName', 'ERA signal')
legend(t2, 'Location', 'southeast', 'Box', 'off')

modeIdx = 1:7;

phaseEra = phaseEra+pi/2;

alpEra = dfEra .* (2*pi.*fEra);

% % --- percentage errors (element-wise, signed) ---
% % guard against divide-by-zero where a true value is 0 (e.g. phase mode 1)
% pctErr = @(est, tru) 100 * (est - tru) ./ tru;
% 
% err_f     = pctErr(fEra(:),      f_true(:));
% err_alpha = pctErr(alpEra(:),    alpha_true(:));
% err_amp   = pctErr(ampEra(:),    amp_true(:));
% err_phase = pctErr(phaseEra(:),  phase_true(:));
% 
% % --- assemble and print ---
% Terr = table(modeIdx(:), f_true(:), err_f, err_alpha, err_amp, err_phase, ...
%     'VariableNames', {'Mode','f_Hz','err_f','err_alpha','err_amp','err_phase'});
% disp(Terr)
% 
% % --- emit LaTeX body rows (tab/ampersand separated, 3 sig figs on errors) ---
% fprintf('\n--- LaTeX rows ---\n');
% for k = 1:numel(modeIdx)
%     fprintf('%d & %g & %.3g & %.3g & %.3g & %.3g \\\\\n', ...
%         modeIdx(k), f_true(k), err_f(k), err_alpha(k), err_amp(k), err_phase(k));
% end

f2 = figure;
til2 = tiledlayout(1, 2);

t1 = nexttile(til2);
hold(t1, "on")
xlabel("Frequency (Hz)")
ylabel("Decay")

t2 = nexttile(til2);
hold(t2, "on")
xlabel('Phase (rad)')
ylabel('Amplitude')

plot(t1, f_true, alpha_true, 'ok', 'LineWidth', 1.5, 'MarkerSize', 10, 'DisplayName', 'Simulated values')
plot(t1, fEra, alpEra, 'rx', 'LineWidth', 1.5,  'MarkerSize', 8, 'DisplayName', 'Simulated values')
legend

plot(t2, phase_true, amp_true, 'ok', 'LineWidth', 1.5, 'MarkerSize', 10, 'DisplayName', 'Simulated values')
plot(t2, phaseEra, ampEra, 'rx', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Simulated values')
legend
%%
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