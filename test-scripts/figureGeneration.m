%% Case 2 - Signals with noise.
clear;
clc;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
msTime = 1000 * t;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

SNR = 5;

freqs = [500, 1000, 1500];

amplitudes = [1.5, 1.0, 0.5];

alpha = 50 : 30 : 110;

phases = 0 : 120 : 240;

M = length(freqs);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amplitudes(m) * exp(-alpha(m) * t) .* ...
        sin(2 * pi * freqs(m) * t + (pi / 180) * phases(m));
end

cleanSig = sum(sigs);

rng(42, 'twister')
%%%%%%%%%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%%%%
Es = sum(cleanSig.^2);
En = Es/(10^(SNR/10));
sign0 = 1 - 2*rand(1,Ns);
En0 = sum(sign0.^2);
an = sqrt(En/En0);
sign = an*sign0;
sig = cleanSig + sign;

opts = eraConfig();
opts.returnDebug = true;
svdTolerance = 3;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, frequencies, decayFactors, modes, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

[~, idx] = sort(frequencies);

frequencies = frequencies(idx);
decayFactors = decayFactors(idx);
modes = modes(idx);
debug.Reconstruction.modalImpulses = debug.Reconstruction.modalImpulses(idx, :);
debug.Modal.decayRates = debug.Modal.decayRates(idx);

f = figure;
til = tiledlayout(2, 2, 'TileSpacing','compact');
setFigureStyle('Style', 'Report')
ylabel(til, 'Amplitude', 'Interpreter', 'latex')
xlabel(til, 'Time (ms)', 'Interpreter', 'latex')

tf = nexttile(1);
title('\textbf{(a)}', 'Interpreter', 'latex')
hold(tf, "on")
ylim([-3 3])
xlim([0 30])
legend

tf2 = nexttile(2);
title('\textbf{(b)}', 'Interpreter', 'latex')
hold(tf2, "on")
ylim([-3 3])
xlim([0 30])
legend

tf3 = nexttile(3);
title('\textbf{(c)}', 'Interpreter', 'latex')
hold(tf3, "on")
ylim([-3 3])
xlim([0 30])
legend

tf4 = nexttile(4);
title('\textbf{(d)}', 'Interpreter', 'latex')
hold(tf4, "on")
ylim([-3 3])
xlim([0 30])
legend

plot(tf, msTime, sig, 'k', 'LineWidth', 1.3, 'DisplayName', 'Noisy impulse')

plot(tf2, msTime, cleanSig, 'LineWidth', 1.5, DisplayName="Clean impulse")
hold(tf2, 'on')
plot(tf2, msTime, outputSignal, '--', 'LineWidth', 1.5, DisplayName="Reconstructed impulse")

plot(tf3, msTime, sigs(1 : end, :), "LineWidth", 1.5)
legend(tf3, {['Frequency: ' num2str(freqs(1))], ...
        ['Frequency: ' num2str(freqs(2))], ...
        ['Frequency: ' num2str(freqs(3))]}, ...
        'Location', 'southeast')

plot(tf4, msTime, debug.Reconstruction.modalImpulses, "LineWidth", 1.5)
frequencies = round(frequencies, 2);
legend(tf4, {['Frequency: ' num2str(frequencies(1))], ...
        ['Frequency: ' num2str(frequencies(2))], ...
        ['Frequency: ' num2str(frequencies(3))]}, ...
        'Location', 'southeast')

f = figure;
setFigureStyle('style', 'Report')
plot(freqs, alpha, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Simulated decay rates')
hold on
plot(frequencies, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay rates')
ylim([40 150])
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Decay rate - $\alpha: (s^{-1}$)', 'Interpreter','latex')
legend('Interpreter','latex', 'Location', 'northwest')

NMSE = sum((cleanSig(:) - ...
    outputSignal(:)) .^ 2) / sum(cleanSig(:) .^ 2);

%%
exportFigure('decay-rates', 'Formats', {'pdf'}, 'OutDir', 'Project/reports/first-report/figures/out')

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
