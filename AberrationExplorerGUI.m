function AberrationExplorerGUI(varargin)
%ABERRATIONEXPLORERGUI  Interactive circular/annular/vignetted-pupil
%   monochromatic-aberration teaching tool (Seidel or Zernike wave-
%   aberration description).
%
%   AberrationExplorerGUI() launches the interactive tool.
%   AberrationExplorerGUI('selftest') runs a battery of numerical
%   self-checks (gradient formulas vs. finite differences, zero-
%   aberration Strehl, pupil-sampling correctness, Zernike-fit
%   consistency, OTF normalization/cutoff) and prints PASS/FAIL to the
%   command window instead of opening the GUI -- useful after any edit
%   to the physics/plotting code below.
%
%   This is an independent, from-scratch implementation built around
%   standard optics theory (Seidel and Zernike wave-aberration
%   polynomials, Fourier-optics PSF/OTF relations via the Wiener-
%   Khinchin theorem, and standard geometric-optics models of the
%   knife-edge/wire/Shack-Hartmann tests -- see e.g. Mahajan, "Optical
%   Imaging and Aberrations", or Goodman, "Introduction to Fourier
%   Optics"). It does not reuse code, layout, or text from any other
%   aberration-visualization program; only the underlying physics
%   (which belongs to the field, not to any one piece of software) is
%   shared.
%
%   FEATURES
%     * 9 classic third-order coefficients (Seidel W-sums or the
%       equivalent low-order Zernike circle polynomials Z0-Z8), each
%       with a slider AND a numeric-entry box.
%     * An optional set of 9 higher-order (4th/5th/6th radial order)
%       terms: trefoil, secondary astigmatism, quadrafoil, secondary
%       coma, and secondary spherical aberration -- toggled on with a
%       checkbox.
%     * Non-circular pupils: annular (central obscuration) or a
%       vignetted/decentered aperture (the intersection of two offset
%       circles), selectable from a menu.
%     * A "field survey" mode that overlays several field heights at
%       once on the field-dependent plots (wavefront lineouts, ray
%       intercept, spot diagram, MTF lineouts).
%     * A simple polychromatic PSF/MTF model: an incoherent sum of the
%       monochromatic PSF over several wavelengths spanning a
%       user-set fractional bandwidth.
%     * All 17 of the classic visualizations: wavefront (3D surface,
%       2D map, x/y lineouts), ray-intercept plot, spot diagram, MTF
%       (2D and lineouts), PSF, single/double-pass interferograms,
%       x/y lateral-shearing interferograms, Shack-Hartmann spots, and
%       knife-edge / wire shadow tests (horizontal and vertical).
%
%   USAGE
%       AberrationExplorerGUI
%
%   Requires base MATLAB only (no toolboxes). Tested for syntax and
%   logic via an equivalent Python/NumPy re-implementation of the
%   numerical core (wavefront, gradient, PSF, and OTF formulas were
%   cross-checked there); it has not been run inside MATLAB itself, so
%   please flag anything that doesn't behave as expected.
%
%   Sign / axis convention: pupil coordinates (x,y) are normalized so
%   the pupil edge sits at radius 1 (x^2+y^2=1); the object/image field
%   point is taken to lie along +y at normalized field height h.
%
%   ACKNOWLEDGMENT / PROVENANCE
%     The concept for a tool like this -- an interactive Seidel/Zernike
%     wavefront explorer covering wavefront maps, ray-intercept and
%     spot-diagram plots, MTF/PSF, interferograms, shearing
%     interferograms, Shack-Hartmann spots, and knife-edge/wire tests --
%     is drawn from a MATLAB teaching tool written by Prof. Brian
%     Kruschwitz for optics coursework at the University of Rochester
%     (Institute of Optics). This file is a completely new, independent
%     implementation: no source code, GUI layout, comments, or text from
%     that tool are reused anywhere here, and it adds a substantial set
%     of features (higher-order aberrations, non-circular pupils, field
%     surveys, a polychromatic PSF model, and slider controls) that the
%     original did not have. The physics itself -- Seidel and Zernike
%     wave-aberration theory and the Fourier-optics PSF/OTF relations --
%     is standard material covered in any optics textbook (e.g. Mahajan,
%     "Optical Imaging and Aberrations"; Goodman, "Introduction to
%     Fourier Optics") and is not anyone's original expression, but full
%     credit for the pedagogical concept and the choice of this
%     particular set of visualizations belongs to Prof. Kruschwitz.

% (no output argument -- calling AberrationExplorerGUI() at the command
% line without a semicolon will not print an "ans = ..." echo)

if nargin >= 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'selftest')
    RunSelfTests();
    return;
end

%% ------------------------------------------------------------------
%  Shared GUI state (visible to every nested callback/plot function)
%% ------------------------------------------------------------------
SeidelFlag      = true;     % true = Seidel coefficients, false = Zernike
ShowHO          = false;    % show/using higher-order coefficient panel
FieldSurveyFlag = false;    % overlay several field heights
PolyChromFlag   = false;    % polychromatic PSF/MTF
PupilShape      = 1;        % 1 circular, 2 annular, 3 vignetted/off-axis
PupilParam      = 0.30;     % obscuration ratio (2) or decenter (3)
hField          = 0;        % normalized field height (Seidel mode only)
vCoef           = zeros(1,18);   % [primary(1:9), higher-order(10:18)]
FieldList       = [0 0.5 0.7 1.0];
Bandwidth       = 0.02;     % fractional dlambda/lambda for polychromatic
PlotType        = 1;
GWmat           = [];       % cached [Ex,Ey] gradient maps (for KE/wire/SH)
AnimatingFlag   = false;    % true while a field-height sweep is playing

PrimaryLabelSeidel  = {'W011x (Tilt X)','W011y (Tilt Y)','W111 (Field Tilt)', ...
                       'W020 (Defocus)','W040 (Spherical)','W131 (Coma)', ...
                       'W222 (Astig.)','W220 (Field Curv.)','W311 (Distortion)'};
PrimaryLabelZernike = {'Z0 (Piston)','Z1 (X-Tilt)','Z2 (Y-Tilt)','Z3 (Defocus)', ...
                       'Z4 (Astig. 0°)', 'Z5 (Astig. 45°)', 'Z6 (Coma X)', ...
                       'Z7 (Coma Y)','Z8 (Spherical)'};
HOLabelStr = {'Trefoil X','Trefoil Y','Sec. Astig. 0°','Sec. Astig. 45°', ...
              'Quadrafoil X','Quadrafoil Y','Sec. Coma X','Sec. Coma Y', ...
              'Sec. Spherical'};

PlotTypeList = {'Wavefront (3D Surface)','Wavefront (2D Color Map)', ...
    'Wavefront (X,Y Lineouts)','Ray Intercept Plot','Spot Diagram', ...
    'MTF (2D Map)','MTF (X,Y Lineouts)','PSF', ...
    'Interferogram (Single-Pass)','Interferogram (Double-Pass)', ...
    'Lateral Shear Interferogram (Y-Shear)','Lateral Shear Interferogram (X-Shear)', ...
    'Shack-Hartmann Spots','Knife-Edge Test (Horizontal Edge)', ...
    'Knife-Edge Test (Vertical Edge)','Wire Test (Horizontal Wire)', ...
    'Wire Test (Vertical Wire)','Zernike Decomposition (Bar Chart)'};

FieldSurveyPlots = [3 4 5 7];   % plot-type indices that honor FieldSurveyFlag

% Preset library (primary Seidel coefficients only; higher-order terms and
% field height default to 0 unless given). Field-dependent terms (coma,
% astigmatism, field curvature, distortion) need h ~= 0 to show any effect,
% so presets that use them also set a representative field height.
PresetNames = {'-- Select a Preset --', 'Diffraction-Limited (Reset)', ...
    'Pure Defocus (1 Wave)', 'Pure Spherical (1 Wave)', ...
    'Pure Coma (1 Wave, h=1)', 'Pure Astigmatism (1 Wave, h=1)', ...
    'Coma + Spherical Mix', 'Strong Spherical (2 Waves)', ...
    'Illustrative "Typical Eye" Mix'};
%                 a1 a2 a3   a4   a5   a6   a7   a8 a9    h
PresetTable = [    0  0  0    0    0    0    0    0  0    0   ;   % Reset
                   0  0  0    1    0    0    0    0  0    0   ;   % Defocus
                   0  0  0    0    1    0    0    0  0    0   ;   % Spherical
                   0  0  0    0    0    1    0    0  0    1   ;   % Coma
                   0  0  0    0    0    0    1    0  0    1   ;   % Astig
                   0  0  0    0  0.5  0.8    0    0  0    1   ;   % Coma+Sph
                   0  0  0    0    2    0    0    0  0    0   ;   % Strong Sph
                   0  0  0 0.25 0.05 0.05 0.15    0  0    1  ];   % "Eye"

%% ------------------------------------------------------------------
%  Build the figure and the two panels
%% ------------------------------------------------------------------
bg          = [0.945 0.945 0.955];   % panel / control background
figBg       = [0.90  0.90  0.91];    % figure background (frames the panels)
accentColor = [0.17  0.44  0.70];    % update button / highlight blue
titleColor  = [0.13  0.15  0.24];    % panel title text

SVfig = figure('Units','normalized','Position',[0.05 0.06 0.90 0.88], ...
    'NumberTitle','off', ...
    'Name','Aberration Explorer  (circular / annular / vignetted pupil)', ...
    'Color',figBg,'Tag','AberrationExplorerFig');

InputPanel = uipanel('Parent',SVfig,'Units','normalized', ...
    'Position',[0.010 0.010 0.460 0.980], 'Title','  Input Panel', ...
    'FontSize',12,'FontWeight','bold','ForegroundColor',titleColor, ...
    'BackgroundColor',bg,'BorderType','line');

OutputPanel = uipanel('Parent',SVfig,'Units','normalized', ...
    'Position',[0.485 0.010 0.505 0.980], 'Title','  Output Panel', ...
    'FontSize',12,'FontWeight','bold','ForegroundColor',titleColor, ...
    'BackgroundColor',bg,'BorderType','line');

%% ---- Input panel: header controls ---------------------------------
uicontrol(InputPanel,'Style','text','String','Coefficient Type:', ...
    'Units','normalized','Position',[0.02 0.947 0.30 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',9);
PolyTypeMenu = uicontrol(InputPanel,'Style','popupmenu', ...
    'String',{'Seidel','Zernike'},'Units','normalized', ...
    'Position',[0.33 0.943 0.28 0.036],'FontSize',9,'Callback',@onPolyTypeChange);

hLabel = uicontrol(InputPanel,'Style','text','String','Field Height h:', ...
    'Units','normalized','Position',[0.63 0.947 0.20 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',9);
hEdit = uicontrol(InputPanel,'Style','edit','String','0', ...
    'Units','normalized','Position',[0.84 0.943 0.13 0.036],'FontSize',9, ...
    'Callback',@onHEdit);

uicontrol(InputPanel,'Style','text','String','Pupil:', ...
    'Units','normalized','Position',[0.02 0.900 0.10 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',9);
PupilMenu = uicontrol(InputPanel,'Style','popupmenu', ...
    'String',{'Circular','Annular','Vignetted (Off-Axis)'}, ...
    'Units','normalized','Position',[0.13 0.896 0.34 0.036],'FontSize',9, ...
    'Callback',@onPupilChange);
PupilParamLabel = uicontrol(InputPanel,'Style','text','String','', ...
    'Units','normalized','Position',[0.49 0.900 0.22 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',8);
PupilParamEdit = uicontrol(InputPanel,'Style','edit','String','0.30', ...
    'Units','normalized','Position',[0.72 0.896 0.13 0.036],'FontSize',9, ...
    'Enable','off','Callback',@onPupilParamEdit);

HOCheck = uicontrol(InputPanel,'Style','checkbox', ...
    'String','Show Higher-Order (4th/5th/6th) Terms', ...
    'Units','normalized','Position',[0.02 0.862 0.94 0.030], ...
    'BackgroundColor',bg,'FontSize',9,'Callback',@onHOToggle);
HONoteText = uicontrol(InputPanel,'Style','text', ...
    'String','Note: shown without additional field (h) dependence -- a simplification; see source comments.', ...
    'Units','normalized','Position',[0.04 0.840 0.94 0.020], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',7, ...
    'FontAngle','italic','ForegroundColor',[0.40 0.40 0.40],'Visible','off');

FieldSurveyCheck = uicontrol(InputPanel,'Style','checkbox', ...
    'String','Field Survey (Seidel Only), h List:', ...
    'Units','normalized','Position',[0.02 0.801 0.46 0.030], ...
    'BackgroundColor',bg,'FontSize',8.5,'Callback',@onFieldSurveyToggle);
FieldListEdit = uicontrol(InputPanel,'Style','edit','String','0, 0.5, 0.7, 1.0', ...
    'Units','normalized','Position',[0.50 0.803 0.47 0.030],'FontSize',8, ...
    'Callback',@onFieldListEdit);

PolyChromCheck = uicontrol(InputPanel,'Style','checkbox', ...
    'String','Polychromatic PSF/MTF, Bandwidth dλ/λ =', ...
    'Units','normalized','Position',[0.02 0.762 0.66 0.030], ...
    'BackgroundColor',bg,'FontSize',8.5,'Callback',@onPolyChromToggle);
BandwidthEdit = uicontrol(InputPanel,'Style','edit','String','0.02', ...
    'Units','normalized','Position',[0.70 0.764 0.14 0.030],'FontSize',8, ...
    'Callback',@onBandwidthEdit);

uicontrol(InputPanel,'Style','text','String','Presets:', ...
    'Units','normalized','Position',[0.02 0.719 0.11 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',9);
PresetMenu = uicontrol(InputPanel,'Style','popupmenu','String',PresetNames, ...
    'Units','normalized','Position',[0.14 0.715 0.58 0.036],'FontSize',8);
uicontrol(InputPanel,'Style','pushbutton','String','Load', ...
    'Units','normalized','Position',[0.74 0.715 0.22 0.036],'FontSize',9, ...
    'Callback',@onLoadPreset);

LiveReadoutText = uicontrol(InputPanel,'Style','text','String','', ...
    'Units','normalized','Position',[0.02 0.680 0.96 0.028], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',8, ...
    'FontWeight','bold','ForegroundColor',titleColor);

%% ---- Input panel: two coefficient columns --------------------------
RowTop = 0.635;
RowBot = 0.110;
NRows  = 9;
RowH   = (RowTop - RowBot) / NRows;

% Within each row: a label strip on top, then a slider+edit strip below
% it separated by a small gap -- sized as explicit fractions of RowH so
% the two strips never overlap regardless of RowH.
LabelFrac = 0.36;
GapFrac   = 0.06;
CtlFrac   = 0.40;

PrimLabelH = gobjects(1,9); PrimSlider = gobjects(1,9); PrimEdit = gobjects(1,9);
HOLabelH   = gobjects(1,9); HOSlider   = gobjects(1,9); HOEdit   = gobjects(1,9);

for k = 1:9
    yTop = RowTop - (k-1)*RowH;
    labelH = LabelFrac*RowH;
    ctlH   = CtlFrac*RowH;
    yLabel = yTop - labelH;
    yCtl   = yLabel - GapFrac*RowH - ctlH;

    PrimLabelH(k) = uicontrol(InputPanel,'Style','text', ...
        'String',PrimaryLabelSeidel{k}, 'Units','normalized', ...
        'Position',[0.02 yLabel 0.44 labelH], ...
        'HorizontalAlignment','left','FontSize',7.5,'BackgroundColor',bg);
    PrimSlider(k) = uicontrol(InputPanel,'Style','slider', ...
        'Min',-10,'Max',10,'Value',0, 'Units','normalized', ...
        'Position',[0.02 yCtl 0.30 ctlH], 'Callback',{@onPrimSlider,k});
    PrimEdit(k) = uicontrol(InputPanel,'Style','edit','String','0', ...
        'Units','normalized','Position',[0.335 yCtl 0.115 ctlH], ...
        'FontSize',8,'Callback',{@onPrimEdit,k});
    try
        addlistener(PrimSlider(k), 'ContinuousValueChange', {@onPrimSlider,k});
    catch
        % Older MATLAB/Octave without this event: falls back to the
        % standard release-triggered Callback above, already wired.
    end

    HOLabelH(k) = uicontrol(InputPanel,'Style','text', ...
        'String',HOLabelStr{k}, 'Units','normalized', ...
        'Position',[0.52 yLabel 0.44 labelH], ...
        'HorizontalAlignment','left','FontSize',7.5,'BackgroundColor',bg, ...
        'Visible','off');
    HOSlider(k) = uicontrol(InputPanel,'Style','slider', ...
        'Min',-10,'Max',10,'Value',0, 'Units','normalized', ...
        'Position',[0.52 yCtl 0.30 ctlH], 'Callback',{@onHOSlider,k}, ...
        'Visible','off');
    HOEdit(k) = uicontrol(InputPanel,'Style','edit','String','0', ...
        'Units','normalized','Position',[0.835 yCtl 0.115 ctlH], ...
        'FontSize',8,'Callback',{@onHOEdit,k},'Visible','off');
    try
        addlistener(HOSlider(k), 'ContinuousValueChange', {@onHOSlider,k});
    catch
    end
end

uicontrol(InputPanel,'Style','pushbutton','String','Update', ...
    'Units','normalized','Position',[0.34 0.050 0.32 0.042], ...
    'FontSize',11,'FontWeight','bold', ...
    'BackgroundColor',accentColor,'ForegroundColor',[1 1 1], ...
    'Callback',@(~,~) redrawPlot());

uicontrol(InputPanel,'Style','text', ...
    'String','Concept: Prof. B. Kruschwitz, University of Rochester -- Independent Implementation', ...
    'Units','normalized','Position',[0.02 0.010 0.96 0.026], ...
    'HorizontalAlignment','center','BackgroundColor',bg,'FontSize',7, ...
    'FontAngle','italic','ForegroundColor',[0.40 0.40 0.40]);

%% ---- Output panel controls ------------------------------------------
uicontrol(OutputPanel,'Style','text','String','Plot Type:', ...
    'Units','normalized','Position',[0.02 0.958 0.14 0.032], ...
    'HorizontalAlignment','left','BackgroundColor',bg,'FontSize',9);
PlotTypeMenu = uicontrol(OutputPanel,'Style','popupmenu','String',PlotTypeList, ...
    'Units','normalized','Position',[0.16 0.953 0.42 0.045],'FontSize',8.5, ...
    'Callback',@onPlotTypeChange);
AnimateButton = uicontrol(OutputPanel,'Style','pushbutton','String','▶ Animate h', ...
    'Units','normalized','Position',[0.62 0.953 0.36 0.045],'FontSize',9, ...
    'Callback',@onAnimatePlay);

AuxText = uicontrol(OutputPanel,'Style','text','String','', ...
    'Units','normalized','Position',[0.60 0.898 0.38 0.030], ...
    'HorizontalAlignment','center','BackgroundColor',bg,'FontSize',8, ...
    'Visible','off');
AuxSlider = uicontrol(OutputPanel,'Style','slider','Units','normalized', ...
    'Position',[0.60 0.868 0.38 0.028],'Visible','off','Callback',@onAuxSlider);
try
    addlistener(AuxSlider, 'ContinuousValueChange', @onAuxSlider);
catch
end

PlotAxes = axes('Parent',OutputPanel,'Units','normalized', ...
    'OuterPosition',[0.06 0.04 0.90 0.80],'FontSize',8);

updatePupilParamUI();
redrawPlot();

%% =====================================================================
%   Nested callback functions (share the workspace defined above)

%% =====================================================================

    function onPolyTypeChange(src,~)
        SeidelFlag = (get(src,'Value')==1);
        if SeidelFlag
            for kk = 1:9, set(PrimLabelH(kk),'String',PrimaryLabelSeidel{kk}); end
            set(hEdit,'Enable','on'); set(hLabel,'Enable','on');
            set(AnimateButton,'Enable','on');
        else
            for kk = 1:9, set(PrimLabelH(kk),'String',PrimaryLabelZernike{kk}); end
            set(hEdit,'Enable','off'); set(hLabel,'Enable','off');
            set(AnimateButton,'Enable','off');
        end
        redrawPlot();
    end

    function onHEdit(src,~)
        val = str2double(get(src,'String'));
        if ~isnan(val), hField = val; end
        redrawPlot();
    end

    function onLoadPreset(~,~)
        idx = get(PresetMenu,'Value');
        if idx <= 1
            return;   % "-- Select a Preset --" placeholder chosen; do nothing
        end
        row = idx - 1;   % PresetTable has no row for the placeholder entry
        SeidelFlag = true;
        set(PolyTypeMenu,'Value',1);
        for kk = 1:9, set(PrimLabelH(kk),'String',PrimaryLabelSeidel{kk}); end
        set(hEdit,'Enable','on'); set(hLabel,'Enable','on'); set(AnimateButton,'Enable','on');

        vCoef = zeros(1,18);
        vCoef(1:9) = PresetTable(row, 1:9);
        hField = PresetTable(row, 10);

        for kk = 1:9
            set(PrimSlider(kk),'Min',min(-10,vCoef(kk)),'Max',max(10,vCoef(kk)),'Value',vCoef(kk));
            set(PrimEdit(kk),'String',sprintf('%.4g',vCoef(kk)));
            set(HOSlider(kk),'Min',-10,'Max',10,'Value',0);
            set(HOEdit(kk),'String','0');
        end
        set(hEdit,'String',sprintf('%.4g',hField));
        redrawPlot();
    end

    function onAnimatePlay(src,~)
        if AnimatingFlag
            AnimatingFlag = false;   % a running loop below will see this and stop
            return;
        end
        if ~SeidelFlag
            return;   % field height has no effect in Zernike mode
        end
        AnimatingFlag = true;
        set(src,'String','■ Stop');
        prevFieldSurvey = FieldSurveyFlag;
        FieldSurveyFlag = false;   % a static multi-h overlay doesn't apply mid-sweep
        try
            hSteps = linspace(0, 1, 26);
            for hh = hSteps
                if ~AnimatingFlag || ~ishghandle(SVfig)
                    break;
                end
                hField = hh;
                set(hEdit,'String',sprintf('%.3g',hField));
                redrawPlot();
                drawnow;
                pause(0.04);
            end
        catch
            % figure was likely closed mid-animation; fall through to cleanup
        end
        AnimatingFlag = false;
        FieldSurveyFlag = prevFieldSurvey;
        if ishghandle(src)
            set(src,'String','▶ Animate h');
        end
    end

    function updatePupilParamUI()
        switch PupilShape
            case 1
                set(PupilParamLabel,'String',''); set(PupilParamEdit,'Enable','off');
            case 2
                set(PupilParamLabel,'String','Obsc. Ratio:'); set(PupilParamEdit,'Enable','on');
            case 3
                set(PupilParamLabel,'String','Decenter:'); set(PupilParamEdit,'Enable','on');
        end
    end

    function onPupilChange(src,~)
        PupilShape = get(src,'Value');
        updatePupilParamUI();
        redrawPlot();
    end

    function onPupilParamEdit(src,~)
        val = str2double(get(src,'String'));
        if ~isnan(val)
            PupilParam = max(0, min(0.95, val));
            set(src,'String',sprintf('%.4g',PupilParam));  % reflect any clamping visibly
        end
        redrawPlot();
    end

    function onHOToggle(src,~)
        ShowHO = logical(get(src,'Value'));
        if ShowHO, vis = 'on'; else, vis = 'off'; end
        set(HOLabelH,'Visible',vis); set(HOSlider,'Visible',vis); set(HOEdit,'Visible',vis);
        set(HONoteText,'Visible',vis);
        redrawPlot();
    end

    function onFieldSurveyToggle(src,~)
        FieldSurveyFlag = logical(get(src,'Value'));
        redrawPlot();
    end

    function onFieldListEdit(src,~)
        vals = str2num(get(src,'String')); %#ok<ST2NM>
        if isnumeric(vals) && ~isempty(vals), FieldList = vals(:).'; end
        redrawPlot();
    end

    function onPolyChromToggle(src,~)
        PolyChromFlag = logical(get(src,'Value'));
        redrawPlot();
    end

    function onBandwidthEdit(src,~)
        val = str2double(get(src,'String'));
        if ~isnan(val) && val >= 0, Bandwidth = val; end
        redrawPlot();
    end

    function onPrimSlider(src,~,idx)
        vCoef(idx) = get(src,'Value');
        set(PrimEdit(idx),'String',sprintf('%.4g',vCoef(idx)));
        redrawPlot();
    end

    function onPrimEdit(src,~,idx)
        val = str2double(get(src,'String'));
        if ~isnan(val)
            vCoef(idx) = val;
            lo = get(PrimSlider(idx),'Min'); hi = get(PrimSlider(idx),'Max');
            if val < lo, set(PrimSlider(idx),'Min',val); end
            if val > hi, set(PrimSlider(idx),'Max',val); end
            set(PrimSlider(idx),'Value',val);
        end
        redrawPlot();
    end

    function onHOSlider(src,~,idx)
        vCoef(9+idx) = get(src,'Value');
        set(HOEdit(idx),'String',sprintf('%.4g',vCoef(9+idx)));
        redrawPlot();
    end

    function onHOEdit(src,~,idx)
        val = str2double(get(src,'String'));
        if ~isnan(val)
            vCoef(9+idx) = val;
            lo = get(HOSlider(idx),'Min'); hi = get(HOSlider(idx),'Max');
            if val < lo, set(HOSlider(idx),'Min',val); end
            if val > hi, set(HOSlider(idx),'Max',val); end
            set(HOSlider(idx),'Value',val);
        end
        redrawPlot();
    end

    function onPlotTypeChange(src,~)
        PlotType = get(src,'Value');
        redrawPlot();
    end

    function onAuxSlider(src,~)
        applyAux(get(src,'Value'));
    end

%% =====================================================================
%   Drawing dispatcher
%% =====================================================================

    function redrawPlot()
        cla(PlotAxes); view(PlotAxes,2);
        set(PlotAxes,'XTickMode','auto','YTickMode','auto', ...
            'XTickLabelMode','auto','YTickLabelMode','auto','XTickLabelRotation',0);
        try, legend(PlotAxes,'off'); catch, end %#ok<CTCH>
        set(AuxSlider,'Visible','off'); set(AuxText,'Visible','off');

        switch PlotType
            case 1,  plotWFsurf();
            case 2,  plotWFcolor();
            case 3,  plotWFlineouts();
            case 4,  plotRayIntercept();
            case 5,  plotSpotDiagram();
            case 6,  plotMTF2D();
            case 7,  plotMTF1D();
            case 8,  plotPSF();
            case 9,  plotInterferogram(1);
            case 10, plotInterferogram(2);
            case 11
                set(AuxSlider,'Min',0,'Max',1,'Value',0.2,'Visible','on');
                set(AuxText,'String','Shear Fraction = 0.200','Visible','on');
                plotShear(false, 0.2);
            case 12
                set(AuxSlider,'Min',0,'Max',1,'Value',0.2,'Visible','on');
                set(AuxText,'String','Shear Fraction = 0.200','Visible','on');
                plotShear(true, 0.2);
            case 13
                GWmat = computeGradientMap();
                set(AuxSlider,'Min',0,'Max',0.20,'Value',0.05,'Visible','on');
                set(AuxText,'String','(f/#)/Rb = 0.050','Visible','on');
                plotShackHartmann(0.05);
            case 14
                GWmat = computeGradientMap();
                setupAuxFromGrad(false);
                plotKnifeEdge(false, 0);
            case 15
                GWmat = computeGradientMap();
                setupAuxFromGrad(true);
                plotKnifeEdge(true, 0);
            case 16
                GWmat = computeGradientMap();
                setupAuxFromGrad(false);
                plotWireTest(false, 0);
            case 17
                GWmat = computeGradientMap();
                setupAuxFromGrad(true);
                plotWireTest(true, 0);
            case 18
                plotZernikeDecomp();
        end
        updateLiveReadout();
    end

    function updateLiveReadout()
        [Xm,Ym] = meshgrid(linspace(-1,1,101));
        Mask = PupilMaskFn(Xm,Ym,PupilShape,PupilParam);
        Wm = WavefrontMap(vCoef, hField, Xm, Ym);
        Wm(~Mask) = NaN;
        [pv,rmsW] = pvRms(Wm);
        strehlApprox = exp(-(2*pi*rmsW)^2);   % Marechal approximation (small-aberration estimate)
        set(LiveReadoutText, 'String', sprintf(...
            'Current Wavefront:  P-V = %.4f waves   RMS = %.4f waves   Strehl (Marechal approx.) ~ %.4f', ...
            pv, rmsW, strehlApprox));
    end

    function setupAuxFromGrad(vertical)
        if vertical, g = GWmat(:,:,1); else, g = GWmat(:,:,2); end
        finiteVals = g(~isnan(g));
        if isempty(finiteVals), gmax = 1; else, gmax = max(1, max(abs(finiteVals))); end
        set(AuxSlider,'Min',-gmax,'Max',gmax,'Value',0,'Visible','on');
        set(AuxText,'String','Position/Threshold = 0.000','Visible','on');
    end

    function applyAux(val)
        switch PlotType
            case {11,12}
                horiz = (PlotType==12);
                set(AuxText,'String',sprintf('Shear Fraction = %.3f',val));
                plotShear(horiz, val);
            case 13
                set(AuxText,'String',sprintf('(f/#)/Rb = %.3f',val));
                plotShackHartmann(val);
            case {14,15}
                vertical = (PlotType==15);
                set(AuxText,'String',sprintf('KE Threshold = %.3f',val));
                plotKnifeEdge(vertical, val);
            case {16,17}
                vertical = (PlotType==17);
                set(AuxText,'String',sprintf('Wire Position = %.3f',val));
                plotWireTest(vertical, val);
        end
    end

%% =====================================================================
%   Wavefront / gradient helpers that need SeidelFlag from the closure
%% =====================================================================

    function Wm = WavefrontMap(v, hh, Xm, Ym)
        if SeidelFlag
            Wm = SeidelWavefront(v(1:9), hh, Xm, Ym);
        else
            Wm = ZernikeWavefront(v(1:9), Xm, Ym);
        end
        Wm = Wm + HigherOrderWavefront(v(10:18), Xm, Ym);
    end

    function [Ex,Ey] = GradientField(v, hh, Xm, Ym)
        if SeidelFlag
            [Ex,Ey] = SeidelGradient(v(1:9), hh, Xm, Ym);
        else
            [Ex,Ey] = ZernikeGradient(v(1:9), Xm, Ym);
        end
        [Hx,Hy] = HigherOrderGradient(v(10:18), Xm, Ym);
        Ex = Ex + Hx; Ey = Ey + Hy;
    end

    function [pv,rmsW] = pvRms(Wm)
        w = Wm(~isnan(Wm));
        if isempty(w), pv = 0; rmsW = 0; return; end
        pv = max(w) - min(w);
        rmsW = sqrt(mean(w.^2) - mean(w).^2);
    end

%% =====================================================================
%   Individual plot routines
%% =====================================================================

    function plotWFsurf()
        % Note on rendering: an earlier version set the Z-data itself to
        % NaN outside the pupil and handed that straight to surf(). Some
        % renderers draw spiky degenerate faces right at a hard
        % NaN/valid-data boundary (most visible with larger aberration
        % amplitudes, where the surface's slope at the pupil edge is
        % steep) -- this is a known surf()-rendering gotcha, not a
        % computation error. For a circular or annular pupil we sidestep
        % it entirely by meshing in POLAR coordinates (rho in
        % [obscuration, 1], theta in [0, 2*pi]): every grid point then
        % lies exactly inside the pupil by construction, so there is no
        % NaN and no boundary discontinuity for any renderer to
        % mishandle. A vignetted (decentered) pupil has no such simple
        % polar description, so that case still masks a Cartesian grid
        % with NaN at a finer resolution.
        thAx = linspace(0, 2*pi, 141);
        if PupilShape == 3
            % Vignetted pupil: the intersection of two unit circles (one
            % centered at the origin, one at (2*dec,0)) is convex, hence
            % star-shaped about its own interior point (dec,0). The
            % boundary radius measured from that point has a closed
            % form (ray/circle intersection against each of the two
            % circles, keeping the nearer one), so the whole region can
            % still be meshed exactly in polar form with no NaN at all.
            dec = max(0, min(0.95, PupilParam));
            s = sqrt(max(0, 1 - dec^2*sin(thAx).^2));
            rMaxTh = s - dec*abs(cos(thAx));
            rFrac = linspace(0, 1, 70).';
            Rm = rFrac * rMaxTh;
            Tm = repmat(thAx, 70, 1);
            Xm = dec + Rm.*cos(Tm);
            Ym = Rm.*sin(Tm);
        else
            if PupilShape == 2
                rho0 = max(0, min(0.95, PupilParam));
            else
                rho0 = 0;
            end
            rhoAx = linspace(rho0, 1, 70);
            [Rm,Tm] = meshgrid(rhoAx, thAx);
            Xm = Rm.*cos(Tm); Ym = Rm.*sin(Tm);
        end
        Wm = WavefrontMap(vCoef, hField, Xm, Ym);
        [pv,rmsW] = pvRms(Wm);
        surf(PlotAxes, Xm, Ym, Wm, 'EdgeColor','none');
        colormap(PlotAxes,'jet'); view(PlotAxes,3); shading(PlotAxes,'interp');
        xlabel(PlotAxes,'\rho_x'); ylabel(PlotAxes,'\rho_y'); zlabel(PlotAxes,'W (waves)');
        title(PlotAxes, sprintf('Wavefront Surface (P-V=%.4f, RMS=%.4f Waves)',pv,rmsW));
    end

    function plotWFcolor()
        [Xm,Ym] = meshgrid(linspace(-1,1,201));
        Mask = PupilMaskFn(Xm,Ym,PupilShape,PupilParam);
        Wm = WavefrontMap(vCoef, hField, Xm, Ym);
        Wm(~Mask) = NaN;
        [pv,rmsW] = pvRms(Wm);
        imagesc(PlotAxes,[-1 1],[-1 1],Wm); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'jet'); try, colorbar(PlotAxes); catch, end %#ok<CTCH>
        xlabel(PlotAxes,'\rho_x'); ylabel(PlotAxes,'\rho_y');
        title(PlotAxes, sprintf('Wavefront Map (P-V=%.4f, RMS=%.4f Waves)',pv,rmsW));
    end

    function plotWFlineouts()
        v = linspace(-1,1,201);
        hold(PlotAxes,'on');
        if FieldSurveyFlag && SeidelFlag
            hs = FieldList; cmap = lines(numel(hs)); leg = {};
            for ii = 1:numel(hs)
                Wx = WavefrontMap(vCoef, hs(ii), v, zeros(size(v)));
                Wy = WavefrontMap(vCoef, hs(ii), zeros(size(v)), v);
                plot(PlotAxes, v, Wx, '-',  'Color',cmap(ii,:));
                plot(PlotAxes, v, Wy, '--', 'Color',cmap(ii,:));
                leg{end+1} = sprintf('W_x, h=%.2g',hs(ii)); %#ok<AGROW>
                leg{end+1} = sprintf('W_y, h=%.2g',hs(ii)); %#ok<AGROW>
            end
            legend(PlotAxes, leg, 'FontSize',7,'Location','bestoutside');
        else
            Wx = WavefrontMap(vCoef, hField, v, zeros(size(v)));
            Wy = WavefrontMap(vCoef, hField, zeros(size(v)), v);
            plot(PlotAxes, v, Wx, 'b', v, Wy, 'r--');
            legend(PlotAxes, {'W_x','W_y'});
        end
        hold(PlotAxes,'off');
        xlabel(PlotAxes,'\rho_x, \rho_y'); ylabel(PlotAxes,'W (waves)');
        title(PlotAxes,'Wavefront Lineouts');
    end

    function plotRayIntercept()
        v = linspace(-1,1,201);
        hold(PlotAxes,'on'); plot(PlotAxes,[-1 1],[0 0],'k-');
        if FieldSurveyFlag && SeidelFlag
            hs = FieldList; cmap = lines(numel(hs)); leg = {};
            for ii = 1:numel(hs)
                [ex,~] = GradientField(vCoef, hs(ii), v, zeros(size(v)));
                [~,ey] = GradientField(vCoef, hs(ii), zeros(size(v)), v);
                plot(PlotAxes, v, ex, '-',  'Color',cmap(ii,:));
                plot(PlotAxes, v, ey, '--', 'Color',cmap(ii,:));
                leg{end+1} = sprintf('\\epsilon_x, h=%.2g',hs(ii)); %#ok<AGROW>
                leg{end+1} = sprintf('\\epsilon_y, h=%.2g',hs(ii)); %#ok<AGROW>
            end
            legend(PlotAxes, leg, 'FontSize',7,'Location','bestoutside');
        else
            [ex,~] = GradientField(vCoef, hField, v, zeros(size(v)));
            [~,ey] = GradientField(vCoef, hField, zeros(size(v)), v);
            plot(PlotAxes, v, ex, 'b', v, ey, 'r--');
            legend(PlotAxes, {'\epsilon_x','\epsilon_y'});
        end
        hold(PlotAxes,'off');
        xlabel(PlotAxes,'\rho_x, \rho_y'); ylabel(PlotAxes,'\epsilon (\lambda/NA)');
        title(PlotAxes,'Ray Intercept Plot');
    end

    function plotSpotDiagram()
        Nspots = 3000;
        hold(PlotAxes,'on');
        maxr = 1;
        if FieldSurveyFlag && SeidelFlag
            hs = FieldList; cmap = lines(numel(hs)); leg = {};
            for ii = 1:numel(hs)
                [px,py] = SamplePupilPoints(PupilShape,PupilParam,Nspots);
                [ex,ey] = GradientField(vCoef, hs(ii), px, py);
                plot(PlotAxes, ex, ey, '.', 'Color',cmap(ii,:), 'MarkerSize',3);
                maxr = max(maxr, max(hypot(ex,ey)));
                leg{end+1} = sprintf('h=%.2g',hs(ii)); %#ok<AGROW>
            end
            legend(PlotAxes, leg, 'Location','bestoutside');
        else
            [px,py] = SamplePupilPoints(PupilShape,PupilParam,Nspots);
            [ex,ey] = GradientField(vCoef, hField, px, py);
            plot(PlotAxes, ex, ey, 'r.', 'MarkerSize',3);
            maxr = max(maxr, max(hypot(ex,ey)));
        end
        hold(PlotAxes,'off'); axis(PlotAxes,'equal');
        lim = 1.25*maxr;
        xlim(PlotAxes,[-lim lim]); ylim(PlotAxes,[-lim lim]);
        xlabel(PlotAxes,'\epsilon_x (\lambda/NA)'); ylabel(PlotAxes,'\epsilon_y (\lambda/NA)');
        title(PlotAxes,'Spot Diagram');
    end

    function plotMTF2D()
        [PSF,~,~,dRho] = CalcPSF(vCoef,hField,SeidelFlag,PupilShape,PupilParam,PolyChromFlag,Bandwidth);
        [OTF,fn] = CalcOTF(PSF,dRho);
        imagesc(PlotAxes, fn, fn, abs(OTF)); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'jet'); try, colorbar(PlotAxes); catch, end %#ok<CTCH>
        hold(PlotAxes,'on'); th = linspace(0,2*pi,200);
        plot(PlotAxes, 2*cos(th), 2*sin(th), 'w-','LineWidth',1.5);
        hold(PlotAxes,'off'); xlim(PlotAxes,[-2.4 2.4]); ylim(PlotAxes,[-2.4 2.4]);
        xlabel(PlotAxes,'f_x / f_c'); ylabel(PlotAxes,'f_y / f_c'); title(PlotAxes,'MTF');
    end

    function plotMTF1D()
        NptsShared = ChoosePupilSampling(vCoef, hField, SeidelFlag);
        [PSF,~,~,dRho] = CalcPSF(vCoef,hField,SeidelFlag,PupilShape,PupilParam,PolyChromFlag,Bandwidth,NptsShared);
        [OTF,fn] = CalcOTF(PSF,dRho);
        [PSFdl,~,~,dRhoDl] = CalcPSF(zeros(1,18),0,SeidelFlag,PupilShape,PupilParam,false,0,NptsShared);
        [OTFdl,~] = CalcOTF(PSFdl,dRhoDl);
        c = floor(numel(fn)/2) + 1;
        fpos = fn(c:end);
        plot(PlotAxes, fpos, abs(OTF(c,c:end)), 'b', ...
                       fpos, abs(OTF(c:end,c)), 'r--', ...
                       fpos, abs(OTFdl(c,c:end)), 'k:','LineWidth',1.5);
        xlim(PlotAxes,[0 2.2]);
        legend(PlotAxes, {'MTF_x','MTF_y','Diffraction Limit'});
        xlabel(PlotAxes,'f / f_c'); ylabel(PlotAxes,'MTF'); title(PlotAxes,'MTF Lineouts');
    end

    function plotPSF()
        [PSF,fx,fy] = CalcPSF(vCoef,hField,SeidelFlag,PupilShape,PupilParam,PolyChromFlag,Bandwidth);
        strehl = max(PSF(:));
        imagesc(PlotAxes, fx, fy, sqrt(PSF)); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'gray'); try, colorbar(PlotAxes); catch, end %#ok<CTCH>
        hold(PlotAxes,'on'); th = linspace(0,2*pi,200);
        plot(PlotAxes, 0.61*cos(th), 0.61*sin(th), 'r-');
        hold(PlotAxes,'off'); xlim(PlotAxes,[-4 4]); ylim(PlotAxes,[-4 4]);
        xlabel(PlotAxes,'x (\lambda/NA)'); ylabel(PlotAxes,'y (\lambda/NA)');
        title(PlotAxes, sprintf('PSF (Sqrt Shown), Strehl = %.4f', strehl));
    end

    function plotInterferogram(npass)
        [Xm,Ym] = meshgrid(linspace(-1,1,201));
        Mask = PupilMaskFn(Xm,Ym,PupilShape,PupilParam);
        Wm = WavefrontMap(vCoef, hField, Xm, Ym);
        I = 0.5*(1 + cos(2*pi*npass*Wm));
        I(~Mask) = NaN;
        imagesc(PlotAxes,[-1 1],[-1 1],I); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'gray');
        hold(PlotAxes,'on'); th = linspace(0,2*pi,200);
        plot(PlotAxes, cos(th), sin(th), 'r-');
        plot(PlotAxes,[-1 1],[0 0],'r-'); plot(PlotAxes,[0 0],[-1 1],'r-');
        hold(PlotAxes,'off');
        xlabel(PlotAxes,'\rho_x'); ylabel(PlotAxes,'\rho_y');
        if npass == 1, ttl = 'Interferogram (Single-Pass)'; else, ttl = 'Interferogram (Double-Pass)'; end
        title(PlotAxes, ttl);
    end

    function plotShear(horizShear, frac)
        halfExt = 1 + frac/2;
        v = linspace(-halfExt, halfExt, 221);
        [Xm,Ym] = meshgrid(v);
        d = frac/2;
        vNoTilt = vCoef; vNoTilt(1) = 0; vNoTilt(2) = 0;
        if horizShear
            X1 = Xm-d; Y1 = Ym; X2 = Xm+d; Y2 = Ym;
            W1 = WavefrontMap(vNoTilt, hField, X1, Y1);
            W2 = WavefrontMap(vNoTilt, hField, X2, Y2) + vCoef(2)*Ym;
        else
            X1 = Xm; Y1 = Ym-d; X2 = Xm; Y2 = Ym+d;
            W1 = WavefrontMap(vNoTilt, hField, X1, Y1);
            W2 = WavefrontMap(vNoTilt, hField, X2, Y2) + vCoef(1)*Xm;
        end
        M1 = PupilMaskFn(X1,Y1,PupilShape,PupilParam);
        M2 = PupilMaskFn(X2,Y2,PupilShape,PupilParam);
        Field = double(M1).*exp(2i*pi*W1) + double(M2).*exp(2i*pi*W2);
        I = abs(Field).^2;
        imagesc(PlotAxes, v, v, I); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'gray');
        xlabel(PlotAxes,'\rho_x'); ylabel(PlotAxes,'\rho_y');
        if horizShear, ttl = 'X-Shearing Interferogram'; else, ttl = 'Y-Shearing Interferogram'; end
        title(PlotAxes, ttl);
    end

    function G = computeGradientMap()
        [Xm,Ym] = meshgrid(linspace(-1,1,201));
        Mask = PupilMaskFn(Xm,Ym,PupilShape,PupilParam);
        [Ex,Ey] = GradientField(vCoef, hField, Xm, Ym);
        Ex(~Mask) = NaN; Ey(~Mask) = NaN;
        G = cat(3, Ex, Ey);
    end

    function plotShackHartmann(fno)
        spacing = 0.4; extent = 0.9;
        coordsAx = -extent:spacing:extent;
        [gx,gy] = meshgrid(coordsAx);
        keep = hypot(gx,gy) <= extent + 1e-9;
        refX = gx(keep).'; refY = gy(keep).';
        Ex = GWmat(:,:,1); Ey = GWmat(:,:,2);
        [Xg,Yg] = meshgrid(linspace(-1,1,size(Ex,1)));
        spotX = zeros(size(refX)); spotY = zeros(size(refY));
        half = spacing/2;
        for ii = 1:numel(refX)
            m = (abs(Xg-refX(ii)) < half) & (abs(Yg-refY(ii)) < half) & ~isnan(Ex);
            if any(m(:))
                mgx = mean(Ex(m)); mgy = mean(Ey(m));
            else
                mgx = 0; mgy = 0;
            end
            spotX(ii) = refX(ii) - fno*mgx;
            spotY(ii) = refY(ii) - fno*mgy;
        end
        plot(PlotAxes, refX, refY, 'k+', spotX, spotY, 'ro');
        hold(PlotAxes,'on'); th = linspace(0,2*pi,200);
        plot(PlotAxes, cos(th), sin(th), 'r-','LineWidth',1.5);
        hold(PlotAxes,'off'); axis(PlotAxes,[-1.1 1.1 -1.1 1.1]); axis(PlotAxes,'square');
        legend(PlotAxes, {'Reference','Test'});
        xlabel(PlotAxes,'\rho_x'); ylabel(PlotAxes,'\rho_y'); title(PlotAxes,'Shack-Hartmann Spots');
    end

    function plotKnifeEdge(vertical, thr)
        if vertical, g = GWmat(:,:,1); else, g = GWmat(:,:,2); end
        valid = ~isnan(g);
        shadow = double(valid);
        shadow(g < thr) = 0;
        shadow = flipud(fliplr(shadow));
        imagesc(PlotAxes,[-1 1],[-1 1],shadow); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'gray');
        hold(PlotAxes,'on'); th = linspace(0,2*pi,200);
        plot(PlotAxes, cos(th), sin(th), 'r--','LineWidth',1.5);
        hold(PlotAxes,'off');
        if vertical, xl = '-\rho_x'; yl = '-\rho_y'; else, xl = '\rho_x'; yl = '\rho_y'; end
        xlabel(PlotAxes,xl); ylabel(PlotAxes,yl);
        title(PlotAxes, sprintf('Knife-Edge Test, Threshold = %.3f', thr));
    end

    function plotWireTest(vertical, pos)
        wireWidth = 0.2;
        if vertical, g = GWmat(:,:,1); else, g = GWmat(:,:,2); end
        valid = ~isnan(g);
        shadow = double(valid);
        band = (g > pos-wireWidth/2) & (g < pos+wireWidth/2);
        shadow(band) = 0;
        shadow = flipud(fliplr(shadow));
        imagesc(PlotAxes,[-1 1],[-1 1],shadow); set(PlotAxes,'YDir','normal');
        axis(PlotAxes,'square'); colormap(PlotAxes,'gray');
        if vertical, xl = '-\rho_x'; yl = '-\rho_y'; else, xl = '\rho_x'; yl = '\rho_y'; end
        xlabel(PlotAxes,xl); ylabel(PlotAxes,yl);
        title(PlotAxes, sprintf('Wire Test, Position = %.3f', pos));
    end

    function plotZernikeDecomp()
        % Least-squares fit of the CURRENT wavefront (whatever mode, field
        % height, and pupil shape are active) onto the 18-term Zernike-
        % style basis (9 classic low-order terms + the 9 higher-order
        % terms), evaluated over the actual pupil shape. Over a plain
        % circular pupil the classic Zernike terms are mutually orthogonal
        % and this recovers the usual per-mode coefficients; over an
        % annular or vignetted pupil the terms are no longer orthogonal,
        % but the least-squares fit is still a well-defined "best
        % representation in this basis over this aperture" and is
        % labeled as such.
        Nfit = 4000;
        [px,py] = SamplePupilPoints(PupilShape, PupilParam, Nfit);
        w = WavefrontMap(vCoef, hField, px, py).';
        A = [ZernikeBasisMatrix(px,py), HigherOrderBasisMatrix(px,py)];
        c = A\w;
        resid = A*c - w;
        resRMS = sqrt(mean(resid.^2));
        wRMS = sqrt(mean((w-mean(w)).^2));
        if wRMS > 1e-12
            pctResid = 100*resRMS/wRMS;
        else
            pctResid = 0;
        end

        allLabels = [PrimaryLabelZernike, HOLabelStr];
        bar(PlotAxes, 1:18, c);
        set(PlotAxes,'XTick',1:18,'XTickLabel',allLabels,'XTickLabelRotation',60);
        grid(PlotAxes,'on');
        ylabel(PlotAxes,'Fitted Coefficient (Waves)');
        shapeNote = '';
        if PupilShape ~= 1
            shapeNote = ' (non-circular pupil: least-squares fit, not an orthogonal mode decomposition)';
        end
        title(PlotAxes, sprintf('Zernike Decomposition -- Fit Residual = %.2f%% of Wavefront RMS%s', ...
            pctResid, shapeNote), 'FontSize',9);
    end

end  % ==================== end of AberrationExplorerGUI ====================


%% =========================================================================
%   Local (non-nested) computational functions
%%   -- standard Seidel / Zernike wave-aberration theory and standard
%%      Fourier-optics diffraction relations; see file header for
%%      references. All partial derivatives below were derived
%%      symbolically (and cross-checked numerically) rather than
%%      transcribed from any other source.
%% =========================================================================

function Wm = SeidelWavefront(a, h, Xm, Ym)
%SEIDELWAVEFRONT  Third-order Seidel wave-aberration sum, in waves.
%   a = [W011x W011y W111 W020 W040 W131 W222 W220 W311]
    rho2 = Xm.^2 + Ym.^2;
    Wm = a(1)*Xm + a(2)*Ym + a(3)*h.*Ym ...
       + a(4)*rho2 + a(5)*rho2.^2 ...
       + a(6)*h.*Ym.*rho2 ...
       + a(7)*h.^2.*Ym.^2 + a(8)*h.^2.*rho2 ...
       + a(9)*h.^3.*Ym;
end

function [Ex,Ey] = SeidelGradient(a, h, Xm, Ym)
%SEIDELGRADIENT  Transverse ray-aberration components eps = -dW/d(coord)
%   for the Seidel wavefront above (normalized by lambda/NA).
    Ex = -( a(1) + 2*a(4)*Xm + 4*a(5)*Xm.*(Xm.^2+Ym.^2) ...
            + 2*a(6)*h.*Xm.*Ym + 2*a(8)*h.^2.*Xm );
    Ey = -( a(2) + a(3)*h + 2*a(4)*Ym + 4*a(5)*Ym.*(Xm.^2+Ym.^2) ...
            + a(6)*h.*(Xm.^2 + 3*Ym.^2) ...
            + 2*a(7)*h.^2.*Ym + 2*a(8)*h.^2.*Ym + a(9)*h.^3 );
end

function Wm = ZernikeWavefront(a, Xm, Ym)
%ZERNIKEWAVEFRONT  Standard low-order (0-8) circle polynomials, in
%   Cartesian pupil coordinates. a = [Z0 Z1 ... Z8] coefficients.
    rho2 = Xm.^2 + Ym.^2;
    Wm = a(1) ...
       + a(2)*Xm + a(3)*Ym ...
       + a(4)*(2*rho2-1) ...
       + a(5)*(Xm.^2-Ym.^2) + a(6)*(2*Xm.*Ym) ...
       + a(7)*(3*rho2-2).*Xm + a(8)*(3*rho2-2).*Ym ...
       + a(9)*(6*rho2.^2 - 6*rho2 + 1);
end

function [Ex,Ey] = ZernikeGradient(a, Xm, Ym)
%ZERNIKEGRADIENT  eps = -dW/d(coord) for the Zernike sum above.
    rho2 = Xm.^2 + Ym.^2;
    dWdx = a(2) + 4*a(4)*Xm + 2*a(5)*Xm + 2*a(6)*Ym ...
         + a(7).*(9*Xm.^2 + 3*Ym.^2 - 2) + 6*a(8)*Xm.*Ym ...
         + a(9).*(24*Xm.^3 + 24*Xm.*Ym.^2 - 12*Xm);
    dWdy = a(3) + 4*a(4)*Ym - 2*a(5)*Ym + 2*a(6)*Xm ...
         + 6*a(7)*Xm.*Ym + a(8).*(3*Xm.^2 + 9*Ym.^2 - 2) ...
         + a(9).*(24*Xm.^2.*Ym + 24*Ym.^3 - 12*Ym);
    Ex = -dWdx; Ey = -dWdy;
end

function Wm = HigherOrderWavefront(v, Xm, Ym)
%HIGHERORDERWAVEFRONT  9 extra 4th/5th/6th-radial-order circle
%   polynomials (trefoil, secondary astigmatism, quadrafoil, secondary
%   coma, secondary spherical), applied as field-independent additions
%   to whichever primary (Seidel or Zernike) wavefront is in use -- a
%   simplification appropriate for a teaching tool, documented here
%   rather than claiming full fifth-order Seidel field theory.
%   v = [trefoil_x trefoil_y secAstig0 secAstig45 quad_x quad_y ...
%        secComa_x secComa_y secSpherical]
    rho2 = Xm.^2 + Ym.^2;
    Wm = v(1)*(Xm.^3 - 3*Xm.*Ym.^2) ...
       + v(2)*(3*Xm.^2.*Ym - Ym.^3) ...
       + v(3)*(4*rho2-3).*(Xm.^2-Ym.^2) ...
       + v(4)*(4*rho2-3).*(2*Xm.*Ym) ...
       + v(5)*(Xm.^4 - 6*Xm.^2.*Ym.^2 + Ym.^4) ...
       + v(6)*(4*Xm.^3.*Ym - 4*Xm.*Ym.^3) ...
       + v(7)*(10*rho2.^2 - 12*rho2 + 3).*Xm ...
       + v(8)*(10*rho2.^2 - 12*rho2 + 3).*Ym ...
       + v(9)*(20*rho2.^3 - 30*rho2.^2 + 12*rho2 - 1);
end

function [Ex,Ey] = HigherOrderGradient(v, Xm, Ym)
%HIGHERORDERGRADIENT  eps = -dW/d(coord) for HigherOrderWavefront above.
    dWdx = v(1)*(3*Xm.^2 - 3*Ym.^2) ...
         + v(2)*(6*Xm.*Ym) ...
         + v(3)*(16*Xm.^3 - 6*Xm) ...
         + v(4)*(24*Xm.^2.*Ym + 8*Ym.^3 - 6*Ym) ...
         + v(5)*(4*Xm.^3 - 12*Xm.*Ym.^2) ...
         + v(6)*(12*Xm.^2.*Ym - 4*Ym.^3) ...
         + v(7)*(50*Xm.^4 + 60*Xm.^2.*Ym.^2 - 36*Xm.^2 + 10*Ym.^4 - 12*Ym.^2 + 3) ...
         + v(8)*(40*Xm.^3.*Ym + 40*Xm.*Ym.^3 - 24*Xm.*Ym) ...
         + v(9)*(120*Xm.^5 + 240*Xm.^3.*Ym.^2 - 120*Xm.^3 + 120*Xm.*Ym.^4 - 120*Xm.*Ym.^2 + 24*Xm);

    dWdy = v(1)*(-6*Xm.*Ym) ...
         + v(2)*(3*Xm.^2 - 3*Ym.^2) ...
         + v(3)*(-16*Ym.^3 + 6*Ym) ...
         + v(4)*(8*Xm.^3 + 24*Xm.*Ym.^2 - 6*Xm) ...
         + v(5)*(-12*Xm.^2.*Ym + 4*Ym.^3) ...
         + v(6)*(4*Xm.^3 - 12*Xm.*Ym.^2) ...
         + v(7)*(40*Xm.^3.*Ym + 40*Xm.*Ym.^3 - 24*Xm.*Ym) ...
         + v(8)*(10*Xm.^4 + 60*Xm.^2.*Ym.^2 - 12*Xm.^2 + 50*Ym.^4 - 36*Ym.^2 + 3) ...
         + v(9)*(120*Xm.^4.*Ym + 240*Xm.^2.*Ym.^3 - 120*Xm.^2.*Ym + 120*Ym.^5 - 120*Ym.^3 + 24*Ym);

    Ex = -dWdx; Ey = -dWdy;
end

function A = ZernikeBasisMatrix(Xm, Ym)
%ZERNIKEBASISMATRIX  9-column design matrix of the individual Z0-Z8
%   circle-polynomial basis shapes (unnormalized), evaluated at the
%   given sample points -- the same shapes ZernikeWavefront sums with
%   coefficients, but returned separately (unweighted) as columns for
%   least-squares fitting/decomposition.
    Xm = Xm(:); Ym = Ym(:);
    rho2 = Xm.^2 + Ym.^2;
    A = [ones(size(Xm)), Xm, Ym, (2*rho2-1), (Xm.^2-Ym.^2), (2*Xm.*Ym), ...
         (3*rho2-2).*Xm, (3*rho2-2).*Ym, (6*rho2.^2-6*rho2+1)];
end

function A = HigherOrderBasisMatrix(Xm, Ym)
%HIGHERORDERBASISMATRIX  9-column design matrix of the individual
%   higher-order shapes, matching the terms HigherOrderWavefront sums
%   with coefficients, returned separately as columns for fitting.
    Xm = Xm(:); Ym = Ym(:);
    rho2 = Xm.^2 + Ym.^2;
    A = [ (Xm.^3 - 3*Xm.*Ym.^2), (3*Xm.^2.*Ym - Ym.^3), ...
          (4*rho2-3).*(Xm.^2-Ym.^2), (4*rho2-3).*(2*Xm.*Ym), ...
          (Xm.^4 - 6*Xm.^2.*Ym.^2 + Ym.^4), (4*Xm.^3.*Ym - 4*Xm.*Ym.^3), ...
          (10*rho2.^2 - 12*rho2 + 3).*Xm, (10*rho2.^2 - 12*rho2 + 3).*Ym, ...
          (20*rho2.^3 - 30*rho2.^2 + 12*rho2 - 1) ];
end

function M = PupilMaskFn(Xm, Ym, shape, param)
%PUPILMASKFN  Boolean pupil-transmission mask.
%   shape: 1 circular, 2 annular (central obscuration = param),
%          3 vignetted (two unit circles offset by 2*param) -- a
%          simple standard model of a vignetted/decentered aperture.
    Rm = hypot(Xm, Ym);
    switch shape
        case 2
            obsc = max(0, min(0.95, param));
            M = (Rm <= 1) & (Rm >= obsc);
        case 3
            dec = max(0, min(0.95, param));
            Rm2 = hypot(Xm - 2*dec, Ym);
            M = (Rm <= 1) & (Rm2 <= 1);
        otherwise
            M = Rm <= 1;
    end
end

function [px,py] = SamplePupilPoints(shape, param, Nspots)
%SAMPLEPUPILPOINTS  Exact (rejection-free) uniform-area sampling of
%   Nspots (x,y) points inside the pupil. Each shape is either directly
%   invertible (circular, annular) or star-shaped about an interior
%   point with a closed-form boundary radius R(theta) (vignetted), so
%   sampling r = R(theta)*sqrt(U), theta ~ U(0,2*pi) gives exactly
%   uniform density over the enclosed area with no rejection and no
%   possibility of returning fewer than Nspots points -- unlike simple
%   accept/reject sampling, which can badly under-deliver for a small
%   or thin pupil (e.g. a heavily obscured annulus or a strongly
%   decentered vignetted aperture).
    theta = 2*pi*rand(1,Nspots);
    switch shape
        case 2  % annular: inner radius obsc, outer radius 1
            obsc = max(0, min(0.95, param));
            rho = sqrt(rand(1,Nspots)*(1-obsc^2) + obsc^2);
            px = rho.*cos(theta); py = rho.*sin(theta);
        case 3  % vignetted: star-shaped about the interior point (dec,0)
            dec = max(0, min(0.95, param));
            s = sqrt(max(0, 1 - dec^2*sin(theta).^2));
            rMaxTh = s - dec*abs(cos(theta));
            rho = rMaxTh.*sqrt(rand(1,Nspots));
            px = dec + rho.*cos(theta); py = rho.*sin(theta);
        otherwise  % circular
            rho = sqrt(rand(1,Nspots));
            px = rho.*cos(theta); py = rho.*sin(theta);
    end
end

function Npts = ChoosePupilSampling(v, h, seidelFlag)
%CHOOSEPUPILSAMPLING  Pick a pupil-grid size fine enough that adjacent
%   pixels' phase does not alias in the FFT-based PSF/OTF, even for the
%   large aberration coefficients the sliders now allow (+/-10 waves).
%   A quick coarse probe estimates the wavefront's peak gradient
%   magnitude, then the grid is sized so consecutive samples differ by
%   at most ~pi/4 radians of phase (comfortably under the pi Nyquist
%   limit). Capped so a pathological input can't make this run away.
    probeN = 65;
    [Xp,Yp] = meshgrid(linspace(-1,1,probeN));
    if seidelFlag
        [ex,ey] = SeidelGradient(v(1:9), h, Xp, Yp);
    else
        [ex,ey] = ZernikeGradient(v(1:9), Xp, Yp);
    end
    [hx,hy] = HigherOrderGradient(v(10:18), Xp, Yp);
    gAll = [abs(ex(:)+hx(:)); abs(ey(:)+hy(:))];
    gAll = gAll(isfinite(gAll));
    if isempty(gAll)
        gmax = 0;
    else
        gmax = max(gAll);
    end
    dRhoNeeded = (pi/4) / max(2*pi*gmax, 1e-6);
    NptsNeeded = ceil(2.2/dRhoNeeded) + 1;
    Npts = max(128, min(400, NptsNeeded));
end

function [PSF, fx, fy, dRho] = CalcPSF(v, h, seidelFlag, shape, param, polyFlag, bw, NptsOverride)
%CALCPSF  Diffraction PSF via Fraunhofer propagation (FFT of the
%   generalized pupil function), normalized to the diffraction-limited
%   peak (Strehl-ratio convention). Optionally sums incoherently over a
%   small wavelength band for a simple polychromatic model: for a fixed
%   physical OPD, the number of waves scales as 1/lambda, so a
%   coefficient set specified at the reference wavelength is divided by
%   the relative wavelength for each color sample.
%   NptsOverride (optional): force a specific pupil-grid size instead of
%   picking one adaptively -- used when a caller needs two PSFs (e.g. an
%   aberrated one and a diffraction-limited reference) to come out on
%   identically-sized grids so they can be compared/indexed together.
    if nargin < 8 || isempty(NptsOverride)
        Npts = ChoosePupilSampling(v, h, seidelFlag);
    else
        Npts = NptsOverride;
    end
    padFactor = max(2, round(2048/Npts));   % keep the padded FFT size bounded
    rhoAx = linspace(-1.1, 1.1, Npts);
    dRho = rhoAx(2) - rhoAx(1);
    [Xm, Ym] = meshgrid(rhoAx, rhoAx);
    Mask = PupilMaskFn(Xm, Ym, shape, param);

    if seidelFlag
        Wm = SeidelWavefront(v(1:9), h, Xm, Ym);
    else
        Wm = ZernikeWavefront(v(1:9), Xm, Ym);
    end
    Wm = Wm + HigherOrderWavefront(v(10:18), Xm, Ym);

    if polyFlag && bw > 0
        lambdaRel = linspace(1 - bw/2, 1 + bw/2, 5);
    else
        lambdaRel = 1;
    end

    Npix = Npts * padFactor;
    PSF = zeros(Npix, Npix);
    for k = 1:numel(lambdaRel)
        Wk = Wm ./ lambdaRel(k);
        P = double(Mask) .* exp(-2i*pi*Wk);
        F = fftshift(fft2(P, Npix, Npix));
        PSF = PSF + abs(F).^2;
    end
    PSF = PSF / numel(lambdaRel);

    Pdl = double(Mask);
    Fdl = fftshift(fft2(Pdl, Npix, Npix));
    PSF = PSF / max(abs(Fdl(:)).^2);

    freqAx = ((-floor(Npix/2)):(ceil(Npix/2)-1)) / (Npix*dRho);
    fx = freqAx; fy = freqAx;
end

function [OTF, fnorm] = CalcOTF(PSF, dRho)
%CALCOTF  OTF via the Wiener-Khinchin route (OTF = FT{PSF}), normalized
%   to unity at zero spatial frequency. Because PSF = |FT{pupil}|^2, a
%   second Fourier transform returns (up to normalization) the
%   autocorrelation of the pupil function; by the standard reciprocal-
%   domain relation for two chained DFTs, the resulting frequency axis
%   has exactly the same sample spacing as the ORIGINAL pupil grid
%   (dRho), so the incoherent cutoff falls at the expected value of 2
%   (twice the unit pupil radius) with no further rescaling needed.
    Npix = size(PSF,1);
    OTF = fftshift(fft2(ifftshift(PSF)));
    c = floor(Npix/2) + 1;
    OTF = OTF / OTF(c,c);
    fnorm = ((-floor(Npix/2)):(ceil(Npix/2)-1)) * dRho;
end

function RunSelfTests()
%RUNSELFTESTS  Numerical self-checks for the physics/plotting core.
%   Invoked via AberrationExplorerGUI('selftest'). Prints PASS/FAIL for
%   each check and a final summary; does not open the GUI.
    nPass = 0; nFail = 0;
    fprintf('=== AberrationExplorerGUI self-test ===\n');

    % ---- 1. Gradient formulas vs. finite differences -------------------
    eps_fd = 1e-6;
    Xp = -0.9 + 1.8*rand(1,25); Yp = -0.9 + 1.8*rand(1,25);

    aS = -1 + 2*rand(1,9); hTest = 0.6;
    [ExS,EyS] = SeidelGradient(aS, hTest, Xp, Yp);
    dWdx_fd = (SeidelWavefront(aS,hTest,Xp+eps_fd,Yp) - SeidelWavefront(aS,hTest,Xp-eps_fd,Yp))/(2*eps_fd);
    dWdy_fd = (SeidelWavefront(aS,hTest,Xp,Yp+eps_fd) - SeidelWavefront(aS,hTest,Xp,Yp-eps_fd))/(2*eps_fd);
    errS = max([abs(ExS+dWdx_fd), abs(EyS+dWdy_fd)]);
    [nPass,nFail] = reportCheck('Seidel gradient matches finite difference', errS < 1e-4, errS, nPass, nFail);

    aZ = -1 + 2*rand(1,9);
    [ExZ,EyZ] = ZernikeGradient(aZ, Xp, Yp);
    dWdx_fd = (ZernikeWavefront(aZ,Xp+eps_fd,Yp) - ZernikeWavefront(aZ,Xp-eps_fd,Yp))/(2*eps_fd);
    dWdy_fd = (ZernikeWavefront(aZ,Xp,Yp+eps_fd) - ZernikeWavefront(aZ,Xp,Yp-eps_fd))/(2*eps_fd);
    errZ = max([abs(ExZ+dWdx_fd), abs(EyZ+dWdy_fd)]);
    [nPass,nFail] = reportCheck('Zernike gradient matches finite difference', errZ < 1e-4, errZ, nPass, nFail);

    aH = -1 + 2*rand(1,9);
    [ExH,EyH] = HigherOrderGradient(aH, Xp, Yp);
    dWdx_fd = (HigherOrderWavefront(aH,Xp+eps_fd,Yp) - HigherOrderWavefront(aH,Xp-eps_fd,Yp))/(2*eps_fd);
    dWdy_fd = (HigherOrderWavefront(aH,Xp,Yp+eps_fd) - HigherOrderWavefront(aH,Xp,Yp-eps_fd))/(2*eps_fd);
    errH = max([abs(ExH+dWdx_fd), abs(EyH+dWdy_fd)]);
    [nPass,nFail] = reportCheck('Higher-order gradient matches finite difference', errH < 1e-3, errH, nPass, nFail);

    % ---- 2. Zero-aberration Strehl = 1 for each pupil shape -------------
    v0 = zeros(1,18);
    for shapeTest = 1:3
        [PSF0,~,~,~] = CalcPSF(v0, 0, true, shapeTest, 0.3, false, 0);
        errStrehl = abs(max(PSF0(:)) - 1);
        [nPass,nFail] = reportCheck(sprintf('Zero-aberration Strehl = 1 (pupil shape %d)', shapeTest), ...
            errStrehl < 1e-6, errStrehl, nPass, nFail);
    end

    % ---- 3. OTF normalization and cutoff --------------------------------
    [PSF0,~,~,dRho0] = CalcPSF(v0, 0, true, 1, 0.3, false, 0);
    [OTF0, fn0] = CalcOTF(PSF0, dRho0);
    c0 = floor(numel(fn0)/2) + 1;
    errOTF0 = abs(OTF0(c0,c0) - 1);
    [nPass,nFail] = reportCheck('OTF(0) = 1', errOTF0 < 1e-9, errOTF0, nPass, nFail);

    lineMTF = abs(OTF0(c0, c0:end));
    idxAbove1pct = find(lineMTF > 0.01, 1, 'last');
    cutoffMeasured = fn0(c0 + idxAbove1pct - 1);
    errCutoff = abs(cutoffMeasured - 2.0);
    [nPass,nFail] = reportCheck('Diffraction-limited MTF cutoff near f/fc = 2', ...
        errCutoff < 0.15, errCutoff, nPass, nFail);

    % ---- 4. Exact analytic pupil sampling: always exactly N points, ----
    %         all satisfying the pupil mask, even for extreme parameters
    testCases = {1,0; 2,0.3; 2,0.9; 3,0.4; 3,0.94};
    for ii = 1:size(testCases,1)
        shp = testCases{ii,1}; prm = testCases{ii,2};
        [px,py] = SamplePupilPoints(shp, prm, 2000);
        okCount = (numel(px) == 2000);
        mask = PupilMaskFn(px, py, shp, prm);
        okInside = all(mask);
        [nPass,nFail] = reportCheck(sprintf('Exact pupil sampling (shape=%d, param=%.2f): count', shp, prm), ...
            okCount, numel(px)-2000, nPass, nFail);
        [nPass,nFail] = reportCheck(sprintf('Exact pupil sampling (shape=%d, param=%.2f): all inside', shp, prm), ...
            okInside, sum(~mask), nPass, nFail);
    end

    % ---- 5. Zernike decomposition of a pure Seidel wavefront is exact --
    %         (every Seidel primary polynomial shape lies in the span of
    %         the classic Z0-Z8 basis; see source comments)
    vSeidelOnly = zeros(1,18);
    vSeidelOnly(1:9) = -1 + 2*rand(1,9);
    hFit = 0.7;
    [pxF,pyF] = SamplePupilPoints(1, 0, 3000);
    wFit = SeidelWavefront(vSeidelOnly(1:9), hFit, pxF, pyF).' + HigherOrderWavefront(vSeidelOnly(10:18),pxF,pyF).';
    Afit = [ZernikeBasisMatrix(pxF,pyF), HigherOrderBasisMatrix(pxF,pyF)];
    cFit = Afit\wFit;
    residFit = Afit*cFit - wFit;
    resRMS = sqrt(mean(residFit.^2));
    wRMS = sqrt(mean((wFit-mean(wFit)).^2));
    relResid = resRMS / max(wRMS, 1e-12);
    [nPass,nFail] = reportCheck('Pure-Seidel wavefront exactly spanned by Zernike+HO basis', ...
        relResid < 1e-6, relResid, nPass, nFail);
    hoCoeffMag = max(abs(cFit(10:18)));
    [nPass,nFail] = reportCheck('...with zero fitted higher-order coefficients', ...
        hoCoeffMag < 1e-6, hoCoeffMag, nPass, nFail);

    fprintf('=== %d passed, %d failed ===\n', nPass, nFail);
    if nFail == 0
        fprintf('ALL TESTS PASSED\n');
    else
        fprintf('%d TEST(S) FAILED -- see above\n', nFail);
    end
end

function [nPass,nFail] = reportCheck(name, isOK, errVal, nPass, nFail)
    if isOK
        fprintf('  [PASS] %s (err=%.3g)\n', name, errVal);
        nPass = nPass + 1;
    else
        fprintf('  [FAIL] %s (err=%.3g)\n', name, errVal);
        nFail = nFail + 1;
    end
end
