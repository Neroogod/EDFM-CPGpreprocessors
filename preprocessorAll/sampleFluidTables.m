function tab = sampleFluidTables(fluid, phases, varargin)
%SAMPLEFLUIDTABLES Turn an MRST fluid's function handles into simulator tables.
%
% SYNOPSIS:
%   tab = sampleFluidTables(fluid, phases)
%   tab = sampleFluidTables(fluid, phases, 'pRange', [50 400]*barsa)
%
% This is the step that makes "export the MRST script and run it" possible.
% An MRST fluid is a set of anonymous functions - b(p), mu(p), kr(s) - while
% Eclipse and CMG both want tabulated properties. This function samples the
% handles on a grid and reports how faithfully a table can represent them,
% so a fluid the table form cannot capture is visible rather than silently
% approximated.
%
% PARAMETERS:
%   fluid  - MRST fluid structure (initSimpleADIFluid, initDeckADIFluid, ...).
%   phases - struct with logical fields .water .oil .gas and, optionally,
%            .disgas .vapoil. Usually built by buildFullDeck from the model.
%
% OPTIONAL ('key', value):
%   'pRange'  - [pMin pMax] in Pa for the PVT sampling. Default is taken
%               from 'pRef' below, widened by a factor of 4 either way.
%   'nP'      - number of PVT samples (default 20).
%   'pRef'    - reference pressure in Pa (default 200 barsa or fluid.pRef).
%   'nSat'    - number of saturation samples per relperm table (default 21).
%   'rsMax'   - maximum solution GOR for a live-oil table. Only used when
%               phases.disgas is true.
%
% RETURNS:
%   tab - struct with
%           .p        [nP x 1] sampling pressures (Pa)
%           .water    .B .mu .b (arrays), .pRef .Bref .muRef .cw .cvw
%           .oil      .B .mu .b (arrays), .rs when live
%           .gas      .B .mu .b (arrays)
%           .sw, .krw, .krow, .pcow          water-oil table
%           .sg, .krg, .krog, .pcog          gas-oil table (3-phase only)
%           .rhoS     [rhoW rhoO rhoG] surface densities (kg/m3)
%           .fit      goodness-of-fit diagnostics, see below
%           .phases   echo of the phase flags actually used
%
% THE FIT DIAGNOSTICS MATTER
%   Eclipse's PVTW and CMG's *CW/*CVW describe water with a single reference
%   point plus one compressibility, i.e. they assume b(p) is log-linear. MRST
%   fluids often are (initSimpleADIFluid uses exp(c*(p-pRef))) but need not
%   be. tab.fit.water.maxRelErr reports the worst relative error of that
%   two-parameter form against the sampled b(p) over the range, and likewise
%   for viscosity. Callers should refuse or warn above a threshold rather
%   than ship a deck whose fluid quietly differs from the MRST one;
%   buildFullDeck does exactly that.
%
% SEE ALSO: buildFullDeck, writeCMGDeck, writeEclipseDeck.

    opt = struct('pRange', [], 'nP', 20, 'pRef', [], 'nSat', 21, 'rsMax', []);
    opt = merge_options(opt, varargin{:});

    if ~isfield(phases, 'disgas'), phases.disgas = false; end
    if ~isfield(phases, 'vapoil'), phases.vapoil = false; end

    % ---- pressure sampling range ---------------------------------------
    pRef = opt.pRef;
    if isempty(pRef)
        if isfield(fluid, 'pRef') && isscalar(fluid.pRef)
            pRef = fluid.pRef;
        else
            pRef = 200*barsa;
        end
    end
    pRange = opt.pRange;
    if isempty(pRange)
        pRange = [pRef/4, pRef*4];
    end
    pRange(1) = max(pRange(1), 1*barsa);   % keep away from vacuum
    assert(pRange(2) > pRange(1), 'sampleFluidTables:pRange', ...
        'pRange must be increasing; got [%g %g].', pRange(1), pRange(2));

    p = linspace(pRange(1), pRange(2), opt.nP)';
    tab = struct();
    tab.p = p;
    tab.pRef = pRef;
    tab.phases = phases;
    tab.fit = struct();

    % ---- surface densities ---------------------------------------------
    rhoS = [NaN NaN NaN];
    if isfield(fluid, 'rhoWS'), rhoS(1) = fluid.rhoWS; end
    if isfield(fluid, 'rhoOS'), rhoS(2) = fluid.rhoOS; end
    if isfield(fluid, 'rhoGS'), rhoS(3) = fluid.rhoGS; end
    tab.rhoS = rhoS;

    % ---- water ----------------------------------------------------------
    if phases.water
        b  = localEval(fluid, 'bW',  p);
        mu = localEval(fluid, 'muW', p);
        tab.water = localPhasePVT(p, b, mu, pRef);
        tab.fit.water = localFitConstComp(p, b, mu, pRef);
    end

    % ---- oil ------------------------------------------------------------
    if phases.oil
        if phases.disgas
            % Live oil: b and mu depend on (p, rs). Saturated table only -
            % the undersaturated branches are emitted by the writers from
            % the saturated end points, which is what both simulators expect
            % when no undersaturated data is available.
            rs = localEval(fluid, 'rsSat', p);
            if isempty(opt.rsMax), rsMax = max(rs); else, rsMax = opt.rsMax; end
            rs = min(rs, rsMax);
            b  = localEval2(fluid, 'bO',  p, rs);
            mu = localEval2(fluid, 'muO', p, rs);
            tab.oil = localPhasePVT(p, b, mu, pRef);
            tab.oil.rs = rs;
            tab.oil.live = true;
        else
            b  = localEval(fluid, 'bO',  p);
            mu = localEval(fluid, 'muO', p);
            tab.oil = localPhasePVT(p, b, mu, pRef);
            tab.oil.rs = zeros(size(p));
            tab.oil.live = false;
        end
        tab.fit.oil = localFitConstComp(p, tab.oil.b, tab.oil.mu, pRef);
    end

    % ---- single-phase water ---------------------------------------------
    % Neither IMEX nor Eclipse has a single-phase option, so a water-only
    % MRST model is exported as oil-water with So = 0 everywhere and the oil
    % given the water's own properties. With no oil in place and none
    % injected the oil equations cannot influence the result, so this is a
    % faithful mapping rather than an invented second phase. The flag lets
    % the writers say so in the deck.
    tab.singlePhaseFiller = false;
    if phases.water && ~phases.oil && ~phases.gas
        tab.oil = tab.water;
        tab.oil.rs = zeros(size(p));
        tab.oil.live = false;
        tab.fit.oil = tab.fit.water;
        tab.singlePhaseFiller = true;
    end

    % ---- gas ------------------------------------------------------------
    if phases.gas
        b  = localEval(fluid, 'bG',  p);
        mu = localEval(fluid, 'muG', p);
        tab.gas = localPhasePVT(p, b, mu, pRef);
        tab.fit.gas = localFitConstComp(p, b, mu, pRef);
    end

    % ---- rock compressibility -------------------------------------------
    % MRST carries this as a pore-volume multiplier; both simulators want a
    % compressibility about a reference pressure.
    if isfield(fluid, 'pvMultR')
        pv = localEval(fluid, 'pvMultR', p);
        [cr, relErr] = localLogSlope(p, pv, pRef);
        tab.rock = struct('cr', cr, 'pRef', pRef, 'pvMult', pv, ...
                          'maxRelErr', relErr);
    else
        tab.rock = struct('cr', 0, 'pRef', pRef, 'pvMult', ones(size(p)), ...
                          'maxRelErr', 0);
    end

    % ---- relative permeability and capillary pressure --------------------
    s = linspace(0, 1, opt.nSat)';

    if phases.water && phases.oil
        tab.sw   = s;
        tab.krw  = localEval(fluid, 'krW', s);
        tab.krow = localKrOW(fluid, 1 - s);   % oil relperm at oil saturation
        tab.pcow = localPc(fluid, 'pcOW', s);
    end

    if phases.gas && phases.oil
        tab.sg   = s;
        tab.krg  = localEval(fluid, 'krG', s);
        tab.krog = localKrOG(fluid, 1 - s);
        tab.pcog = localPc(fluid, 'pcOG', s);
    end

    if phases.water && phases.gas && ~phases.oil
        tab.sw  = s;
        tab.krw = localEval(fluid, 'krW', s);
        tab.krg = localEval(fluid, 'krG', 1 - s);
        tab.pcow = localPc(fluid, 'pcOW', s);
    end

    if tab.singlePhaseFiller
        % The filler oil needs a relperm curve for the deck to be well
        % formed. It is never evaluated at So > 0, so a straight line is
        % both sufficient and honest about carrying no information.
        tab.sw   = s;
        tab.krw  = localEval(fluid, 'krW', s);
        tab.krow = 1 - s;
        tab.pcow = localPc(fluid, 'pcOW', s);
    end
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function v = localEval(fluid, name, x)
% Evaluate a fluid function handle, with a clear error if it is missing.
    assert(isfield(fluid, name), 'sampleFluidTables:missing', ...
        ['fluid.%s is required for the requested phase set but is not ', ...
         'present. Supply it, or exclude that phase.'], name);
    f = fluid.(name);
    if isa(f, 'function_handle')
        v = f(x);
    else
        v = repmat(f, size(x));   % constant property
    end
    v = double(value(v));
    v = reshape(v, [], 1);
end


function v = localEval2(fluid, name, x, y)
% Two-argument form, for live-oil b(p, rs) / mu(p, rs).
    assert(isfield(fluid, name), 'sampleFluidTables:missing', ...
        'fluid.%s is required for a live-oil model.', name);
    f = fluid.(name);
    try
        v = f(x, y, true(size(x)));
    catch
        v = f(x, y);
    end
    v = reshape(double(value(v)), [], 1);
end


function s = localPhasePVT(p, b, mu, pRef)
% Package one phase's sampled properties. MRST's b is the INVERSE formation
% volume factor; both simulators tabulate B, so convert once here and keep
% b alongside for the compressibility fits.
    assert(all(b > 0), 'sampleFluidTables:badB', ...
        'Inverse formation volume factor is non-positive over the sampling range.');
    s = struct();
    s.b  = b;
    s.B  = 1 ./ b;
    s.mu = mu;
    s.pRef  = pRef;
    s.Bref  = interp1(p, s.B,  pRef, 'linear', 'extrap');
    s.muRef = interp1(p, s.mu, pRef, 'linear', 'extrap');
end


function f = localFitConstComp(p, b, mu, pRef)
% How well does the two-parameter (reference + constant compressibility)
% form that PVTW and CMG's *CW/*CVW assume reproduce the sampled data?
%
% c   = d(ln b)/dp, the isothermal compressibility. MRST's b is the inverse
%       formation volume factor, so this is exactly the compressibility both
%       Eclipse (PVTW, ROCK) and CMG (*CW, *CO) expect.
% cmu = d(ln mu)/dp, the viscosity pressure coefficient. CMG's *CVW/*CVO and
%       Eclipse's PVTW viscosibility use this same sign convention, and it is
%       zero for the constant-viscosity fluids MRST most often supplies.
    [c, bErr]   = localLogSlope(p, b, pRef);
    [cmu, mErr] = localLogSlope(p, mu, pRef);
    f = struct('c', c, 'cmu', cmu, ...
               'maxRelErrB', bErr, 'maxRelErrMu', mErr, ...
               'maxRelErr', max(bErr, mErr));
end


function [c, maxRelErr] = localLogSlope(p, y, pRef)
% Least-squares slope of log(y) against p, and the worst relative error of
% the resulting exponential model over the sampled range.
    if all(abs(y - y(1)) <= 1e-14*max(abs(y(1)), 1))
        c = 0; maxRelErr = 0; return          % exactly constant
    end
    if any(y <= 0)
        c = 0; maxRelErr = inf; return        % cannot be modelled this way
    end
    A = [ones(numel(p),1), p(:) - pRef];
    coef = A \ log(y(:));
    c = coef(2);
    yFit = exp(coef(1) + c*(p(:) - pRef));
    maxRelErr = max(abs(yFit - y(:)) ./ max(abs(y(:)), eps));
end


function kr = localKrOW(fluid, so)
% Oil relperm in the water-oil system, as a function of oil saturation.
    if isfield(fluid, 'krOW')
        kr = localEval(fluid, 'krOW', so);
    else
        kr = localEval(fluid, 'krO', so);
    end
end


function kr = localKrOG(fluid, so)
    if isfield(fluid, 'krOG')
        kr = localEval(fluid, 'krOG', so);
    else
        kr = localEval(fluid, 'krO', so);
    end
end


function pc = localPc(fluid, name, s)
% Capillary pressure, zero when the fluid does not define it.
    if isfield(fluid, name)
        pc = localEval(fluid, name, s);
    else
        pc = zeros(numel(s), 1);
    end
end
