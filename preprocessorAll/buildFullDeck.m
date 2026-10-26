function deck = buildFullDeck(model, state0, schedule, varargin)
%BUILDFULLDECK Assemble a complete, simulator-neutral deck from an MRST case.
%
% SYNOPSIS:
%   deck = buildFullDeck(model, state0, schedule)
%   deck = buildFullDeck(model, state0, schedule, 'grdecl', grdecl)
%
% The `modules/preprocessor` writers export the reservoir description -
% geometry, rock, fractures, connections - and stop there, deliberately
% inventing no fluid data. This module goes the rest of the way: it takes
% the same triple an MRST run takes (model, state0, schedule) and produces
% everything a commercial simulator needs to execute the case, so the
% exported deck runs directly instead of being finished by hand in Builder.
%
% PARAMETERS:
%   model    - an MRST AD model (GenericBlackOilModel and its relatives).
%              model.G, model.rock and model.fluid are read from it.
%   state0   - initial state: .pressure and .s (and .rs when live-oil).
%   schedule - MRST schedule with .step.val, .step.control and .control(i).W.
%              Pass [] for a static export with no recurrent section.
%
% OPTIONAL ('key', value):
%   'grdecl'    - corner-point description, required for processGRDECL grids.
%   'datum'     - depth [m] of the top of the grid, overriding the grid's z.
%   'simulator' - 'IMEX' (default), 'GEM' or 'STARS'; recorded in metadata
%                 and used by writeCMGDeck.
%   'fractureModel' - 'auto' (default) | 'matrix' | 'EDFM' | 'pEDFM' | 'CPEDFM'
%   'startDate' - [Y M D] for the recurrent section (default [2026 1 1]).
%   'nP','nSat' - PVT / relperm sampling resolution (see sampleFluidTables).
%   'pRange'    - PVT sampling range [pMin pMax] in Pa. Default spans the
%                 initial pressure field widened by a factor of 4, so the
%                 tables always bracket the states the run will visit.
%   'fitTol'    - largest acceptable relative error when a sampled property
%                 is collapsed to the two-parameter form the deck formats
%                 assume (default 1e-3). Exceeded, buildFullDeck warns and
%                 records it in deck.warnings rather than failing silently.
%   'tol'       - geometric tolerance passed through to buildCMGModel.
%
% RETURNS:
%   deck - simulator-neutral struct:
%            .res      reservoir description from buildCMGModel (grid, rock,
%                      fractures, connections, multipliers, wells)
%            .tab      sampled PVT / rock-fluid tables (sampleFluidTables)
%            .init     initialization: pressures, saturations, datum
%            .wells    per-control well definitions with explicit WI
%            .steps    report steps in days and their control index
%            .meta     phases, model class, units, provenance
%            .warnings cell array of fidelity warnings; empty is the good case
%
% WHAT IS AND IS NOT CARRIED OVER
%   Carried: grid geometry, porosity, permeability, net-to-gross, fracture
%   cells and NNCs, transmissibility multipliers, PVT, relative permeability,
%   capillary pressure, rock compressibility, initial pressure and
%   saturations, wells with MRST's own well indices, controls and limits,
%   and the report-step schedule.
%
%   Not carried: anything with no deck equivalent - solver tuning, MRST-side
%   discretization choices, non-tabulable property functions. Where a
%   property must be approximated to fit a deck format, the residual is
%   measured and reported in .warnings; nothing is quietly rounded off.
%
% SEE ALSO: sampleFluidTables, writeCMGDeck, writeEclipseDeck.

    opt = struct('grdecl', [], 'datum', [], 'simulator', 'IMEX', ...
                 'fractureModel', 'auto', 'startDate', [2026 1 1], ...
                 'nP', 20, 'nSat', 21, 'pRange', [], 'fitTol', 1e-3, ...
                 'tol', 1e-5, 'cmgVersion', '2024.30');
    opt = merge_options(opt, varargin{:});

    assert(isprop(model, 'G') || isfield(model, 'G'), 'buildFullDeck:noG', ...
        'model.G is required.');
    G = model.G;
    rock = model.rock;
    fluid = model.fluid;
    assert(~isempty(fluid), 'buildFullDeck:noFluid', ...
        'model.fluid is empty; there is nothing to export as a fluid model.');

    deck = struct();
    deck.warnings = {};

    % ---- phases ---------------------------------------------------------
    phases = struct('water', localFlag(model, 'water'), ...
                    'oil',   localFlag(model, 'oil'), ...
                    'gas',   localFlag(model, 'gas'), ...
                    'disgas', localFlag(model, 'disgas'), ...
                    'vapoil', localFlag(model, 'vapoil'));
    nPh = phases.water + phases.oil + phases.gas;
    assert(nPh >= 1, 'buildFullDeck:noPhases', 'The model declares no phases.');
    assert(~phases.vapoil, 'buildFullDeck:vapoil', ...
        ['Vaporized oil (vapoil) has no table written by this module yet. ', ...
         'Export with vapoil disabled, or extend sampleFluidTables.']);

    % ---- reservoir description (reuse the existing preprocessor) --------
    deck.res = buildCMGModel(G, rock, ...
        'simulator', opt.simulator, 'cmgVersion', opt.cmgVersion, ...
        'model', opt.fractureModel, 'grdecl', opt.grdecl, ...
        'datum', opt.datum, 'tol', opt.tol, 'wells', []);

    % ---- PVT and rock-fluid tables --------------------------------------
    pRange = opt.pRange;
    if isempty(pRange) && ~isempty(state0) && isfield(state0, 'pressure')
        p0 = double(state0.pressure(:));
        p0 = p0(isfinite(p0) & p0 > 0);
        if ~isempty(p0)
            pRange = [min(p0)/4, max(p0)*4];
        end
    end
    pRefGuess = [];
    if ~isempty(state0) && isfield(state0, 'pressure')
        pRefGuess = mean(double(state0.pressure(:)), 'omitnan');
    end

    deck.tab = sampleFluidTables(fluid, phases, ...
        'pRange', pRange, 'nP', opt.nP, 'nSat', opt.nSat, 'pRef', pRefGuess);

    deck.warnings = [deck.warnings, localFitWarnings(deck.tab, opt.fitTol)];

    % ---- initialization --------------------------------------------------
    deck.init = localBuildInit(G, state0, phases, deck.res);

    % ---- wells and schedule ---------------------------------------------
    [deck.wells, deck.steps, schedWarn] = localBuildSchedule(G, ...
        schedule, phases, deck.res);
    deck.warnings = [deck.warnings, schedWarn];

    % ---- metadata --------------------------------------------------------
    deck.meta = struct( ...
        'schema', 'MRST-full-deck-v1', ...
        'modelClass', class(model), ...
        'simulator', deck.res.meta.simulator, ...
        'cmgVersion', deck.res.meta.cmgVersion, ...
        'gridFamily', deck.res.meta.gridFamily, ...
        'fractureModel', deck.res.meta.fractureModel, ...
        'phases', phases, ...
        'nPhases', nPh, ...
        'unitSystem', 'SI (METRIC)', ...
        'startDate', opt.startDate, ...
        'nCells', deck.res.grid.dims, ...
        'nSteps', numel(deck.steps.dtDays), ...
        'nWellControls', numel(deck.wells), ...
        'gravity', localGravityOn(model), ...
        'created', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function v = localFlag(model, name)
    v = false;
    try
        if isprop(model, name) || isfield(model, name)
            v = logical(model.(name));
        end
    catch
        v = false;
    end
end


function g = localGravityOn(model)
    g = true;
    try
        if isprop(model, 'gravity') || isfield(model, 'gravity')
            g = norm(model.gravity) > 0;
        end
    catch
        g = true;
    end
end


function warns = localFitWarnings(tab, fitTol)
% Report every property whose deck representation departs from the MRST
% function by more than fitTol. These are the places the exported run can
% differ from the MRST run, so they are surfaced rather than buried.
    warns = {};
    names = {'water', 'oil', 'gas'};
    for k = 1:numel(names)
        if ~isfield(tab.fit, names{k}), continue; end
        f = tab.fit.(names{k});
        if ~isfinite(f.maxRelErr)
            warns{end+1} = sprintf(['%s: property is not representable by a ', ...
                'constant-compressibility fit (non-positive values sampled). ', ...
                'The tabulated form is used instead; check the deck.'], names{k}); %#ok<AGROW>
        elseif f.maxRelErr > fitTol
            warns{end+1} = sprintf(['%s: constant-compressibility fit differs ', ...
                'from the MRST function by up to %.2e (B %.2e, mu %.2e) over the ', ...
                'sampled pressure range. Formats that take a single ', ...
                'compressibility (Eclipse PVTW, CMG *CW/*CVW) will carry that ', ...
                'error; the tabulated formats (PVDO/PVDG, CMG *PVT) do not.'], ...
                names{k}, f.maxRelErr, f.maxRelErrB, f.maxRelErrMu); %#ok<AGROW>
        end
    end
    if isfield(tab, 'rock') && tab.rock.maxRelErr > fitTol
        warns{end+1} = sprintf(['rock: pore-volume multiplier is not ', ...
            'log-linear in pressure (max rel. error %.2e). A single ROCK / ', ...
            '*CPOR compressibility cannot reproduce it exactly.'], ...
            tab.rock.maxRelErr);
    end
end


function init = localBuildInit(G, state0, phases, res)
% Initialization straight from state0. Explicit cell-by-cell arrays are used
% rather than an equilibration region: state0 is what MRST actually started
% from, and re-deriving it from contacts would introduce a difference the
% user never asked for.
    init = struct();
    nMat = res.grid.nCells;
    if isempty(state0)
        init.mode = 'none';
        init.pressure = [];
        init.s = [];
        return
    end

    p = double(state0.pressure(:));
    init.mode = 'explicit';
    init.pressure = p;
    init.pMean = mean(p, 'omitnan');
    init.pMin = min(p); init.pMax = max(p);

    if isfield(state0, 's') && ~isempty(state0.s)
        init.s = double(state0.s);
    else
        init.s = [];
    end
    init.phases = phases;

    if isfield(state0, 'rs') && phases.disgas
        init.rs = double(state0.rs(:));
    end

    % Datum for the deck's reference depth, and the pressure there. Taken
    % from the shallowest cell of the matrix grid so that a gravity-
    % initialized state exports with a datum pressure consistent with the
    % depths written into the deck.
    init.datum = res.grid.ztop;
    init.zbot  = res.grid.zbot;
    init.pDatum = init.pMean;

    z = G.cells.centroids(:, 3);
    nUse = min([numel(z), numel(p), nMat]);
    if nUse > 0
        [zMin, ix] = min(z(1:nUse));
        init.pDatum = p(ix);
        init.datum  = zMin;
    end
    init.nMatrixCells = nMat;
end


function [wells, steps, warns] = localBuildSchedule(G, schedule, phases, res)
% Flatten an MRST schedule into (a) one well definition per control period
% and (b) a list of report steps. Well indices come from MRST, never from a
% radial model applied to exported geometry.
    warns = {};
    wells = struct('control', {}, 'W', {});
    steps = struct('dtDays', [], 'control', [], 'tEndDays', []);

    if isempty(schedule)
        return
    end

    assert(isfield(schedule, 'step') && isfield(schedule, 'control'), ...
        'buildFullDeck:schedule', ...
        'schedule must have .step and .control (see simpleSchedule).');

    dt = double(schedule.step.val(:));
    ctrl = double(schedule.step.control(:));
    steps.dtDays = convertTo(dt, day);
    steps.control = ctrl;
    steps.tEndDays = cumsum(steps.dtDays);

    nCtrl = numel(schedule.control);
    for c = 1:nCtrl
        W = schedule.control(c).W;
        if isempty(W)
            wells(c).control = c; wells(c).W = struct([]); continue
        end
        wl = repmat(localEmptyWell(), numel(W), 1);
        for k = 1:numel(W)
            wl(k) = localOneWell(W(k), phases, res, G);
        end
        wells(c).control = c;
        wells(c).W = wl;
    end

    % Completions outside the matrix grid cannot be addressed by an (i,j,k)
    % well card on a padded fracture export.
    nMat = res.grid.nCells;
    for c = 1:numel(wells)
        for k = 1:numel(wells(c).W)
            if any(wells(c).W(k).cells > nMat)
                warns{end+1} = sprintf(['well %s completes in cell(s) beyond ', ...
                    'the matrix grid (fracture cells). Those completions are ', ...
                    'exported at their padded addresses; verify them in the ', ...
                    'target simulator.'], wells(c).W(k).name); %#ok<AGROW>
            end
        end
    end
end


function w = localEmptyWell()
    w = struct('name', '', 'type', '', 'val', 0, 'sign', 0, ...
               'isInjector', false, 'cells', [], 'wi', [], 'uba', [], ...
               'compi', [], 'phase', '', 'refDepth', 0, 'radius', 0, ...
               'status', true, 'lims', []);
end


function w = localOneWell(Wk, phases, res, G)
    w = localEmptyWell();
    w.name = char(Wk.name);
    w.type = char(Wk.type);
    w.val  = double(Wk.val);
    w.sign = double(Wk.sign);
    w.isInjector = w.sign > 0;
    w.cells = double(Wk.cells(:));
    w.wi    = double(Wk.WI(:));
    w.refDepth = double(Wk.refDepth);
    if isfield(Wk, 'r'), w.radius = double(Wk.r(1)); end
    if isfield(Wk, 'status') && ~isempty(Wk.status), w.status = logical(Wk.status(1)); end
    if isfield(Wk, 'lims'), w.lims = Wk.lims; end

    % Injected phase, from compi.
    if isfield(Wk, 'compi') && ~isempty(Wk.compi)
        w.compi = double(Wk.compi(:))';
        w.phase = localPhaseFromCompi(w.compi, phases);
    else
        w.phase = 'WATER';
    end

    % Logical addresses for the well cards.
    w.uba = localCellsToUBA(G, res, w.cells);
end


function name = localPhaseFromCompi(compi, phases)
% compi is ordered by the model's active phases.
    names = {};
    if phases.water, names{end+1} = 'WATER'; end
    if phases.oil,   names{end+1} = 'OIL';   end
    if phases.gas,   names{end+1} = 'GAS';   end
    n = min(numel(names), numel(compi));
    [~, ix] = max(compi(1:n));
    name = names{ix};
end


function uba = localCellsToUBA(G, res, cells)
% Map global cell indices to the (i,j,k) addresses of the exported grid.
%
% Structured and corner-point exports keep MRST's logical numbering, so the
% address is ind2sub of the cell's entry in G.cells.indexMap (which accounts
% for inactive cells). The unstructured export is an Nc x 1 x 1 pseudo-grid,
% where the cell index is the address.
    cells = cells(:);
    if strcmp(res.grid.family, 'unstructured')
        uba = [cells, ones(numel(cells), 1), ones(numel(cells), 1)];
        return
    end

    assert(isfield(G, 'cartDims') && ~isempty(G.cartDims), ...
        'buildFullDeck:noCartDims', ...
        ['Grid family "%s" was detected but G has no cartDims, so well ', ...
         'completions cannot be given logical addresses.'], res.grid.family);

    if isfield(G.cells, 'indexMap') && ~isempty(G.cells.indexMap)
        lin = double(G.cells.indexMap(cells));
    else
        lin = cells;
    end
    [i, j, k] = ind2sub(reshape(double(G.cartDims), 1, []), lin);
    uba = [i(:), j(:), k(:)];
end
