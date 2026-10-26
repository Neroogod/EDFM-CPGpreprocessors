function out = readCMGSr3(sr3File, varargin)
%READCMGSR3 Read an IMEX/GEM/STARS SR3 result file back into MATLAB.
%
% SYNOPSIS:
%   out = readCMGSr3('case.sr3')
%
% CMG's SR3 is an HDF5 container, so the results of a deck exported by this
% module can be read straight back and compared against the MRST run that
% produced it. That comparison is the only way to show an export is faithful
% rather than merely syntactically acceptable.
%
% PARAMETERS:
%   sr3File - path to the .sr3 written by the simulator.
%
% OPTIONAL ('key', value):
%   'props' - cell array of spatial property names to read
%             (default {'PRES','SW','SO','SG'}; missing ones are skipped).
%
% RETURNS:
%   out - struct with
%           .timeDays   report times of the spatial snapshots [nSnap x 1]
%           .props      struct of [nCells x nSnap] arrays, one field per
%                       property actually present (PRES in kPa, saturations
%                       as fractions)
%           .wells      struct with .names, .vars, .timeDays and .data
%                       [nWell x nVar x nStep], as written by the simulator
%           .file       the file that was read
%
% Well variable names of interest: BHP, OILRATSC, WATRATSC, LIQRATSC (rates
% at surface conditions, m3/day in an SI deck) and the cumulative OILVOLSC /
% WATVOLSC.
%
% SEE ALSO: writeCMGDeck, buildFullDeck.

    opt = struct('props', {{'PRES', 'SW', 'SO', 'SG'}});
    opt = merge_options(opt, varargin{:});

    assert(isfile(sr3File), 'readCMGSr3:missing', ...
        'SR3 file not found: %s', sr3File);

    out = struct('file', sr3File, 'timeDays', [], 'props', struct(), ...
                 'wells', struct());

    info = h5info(sr3File);
    names = {info.Groups.Name};

    % ---- spatial snapshots ---------------------------------------------
    if any(strcmp(names, '/SpatialProperties'))
        sp = h5info(sr3File, '/SpatialProperties');
        snapNames = {};
        for k = 1:numel(sp.Groups)
            [~, leaf] = fileparts(sp.Groups(k).Name);
            snapNames{end+1} = leaf; %#ok<AGROW>
        end
        snapNames = sort(snapNames);

        % Snapshot index -> report time, via the master time table.
        idxToTime = containers.Map('KeyType', 'double', 'ValueType', 'double');
        try
            mtt = h5read(sr3File, '/General/MasterTimeTable');
            for k = 1:numel(mtt.Index)
                idxToTime(double(mtt.Index(k))) = double(mtt.OffsetInDays(k));
            end
        catch
            % No master table: fall back to the snapshot ordinal below.
        end

        tVec = zeros(numel(snapNames), 1);
        store = struct();
        for s = 1:numel(snapNames)
            gname = ['/SpatialProperties/' snapNames{s}];
            g = h5info(sr3File, gname);
            idx = str2double(snapNames{s});
            if isKey(idxToTime, idx), tVec(s) = idxToTime(idx); else, tVec(s) = s-1; end
            for d = 1:numel(g.Datasets)
                nm = g.Datasets(d).Name;
                if ~any(strcmpi(nm, opt.props)), continue; end
                v = double(h5read(sr3File, [gname '/' nm]));
                v = v(:);
                if ~isfield(store, nm)
                    store.(nm) = nan(numel(v), numel(snapNames));
                end
                n = min(numel(v), size(store.(nm), 1));
                store.(nm)(1:n, s) = v(1:n);
            end
        end
        out.timeDays = tVec;
        out.props = store;
    end

    % ---- well time series ------------------------------------------------
    if any(strcmp(names, '/TimeSeries'))
        try
            % CMG pads these fixed-width HDF5 strings out to their field
            % length. strtrim on the raw h5read result does not always clear
            % it, so normalize through the string type, which strips
            % reliably, and hand back plain cellstr.
            wn = localCleanNames(h5read(sr3File, '/TimeSeries/WELLS/Origins'));
            vn = localCleanNames(h5read(sr3File, '/TimeSeries/WELLS/Variables'));
            ts = double(h5read(sr3File, '/TimeSeries/WELLS/Timesteps'));
            dt = double(h5read(sr3File, '/TimeSeries/WELLS/Data'));
            out.wells = struct( ...
                'names', {wn}, ...
                'vars',  {vn}, ...
                'step',  ts(:), ...
                'data',  dt);
            % Map the well timestep indices onto report times when possible.
            try
                mtt = h5read(sr3File, '/General/MasterTimeTable');
                m = containers.Map(double(mtt.Index), double(mtt.OffsetInDays));
                tw = zeros(numel(ts), 1);
                for k = 1:numel(ts)
                    if isKey(m, ts(k)), tw(k) = m(ts(k)); else, tw(k) = NaN; end
                end
                out.wells.timeDays = tw;
            catch
                out.wells.timeDays = ts(:);
            end
        catch me
            warning('readCMGSr3:wells', ...
                'Well time series could not be read: %s', me.message);
        end
    end
end


function c = localCleanNames(raw)
% Normalize an HDF5 fixed-width string column to a trimmed cellstr.
% Trailing padding may be spaces or NULs depending on how the field was
% written, so strip both explicitly rather than relying on strtrim alone.
    s = string(raw);
    s = s(:);
    s = erase(s, char(0));
    s = strip(s);
    c = cellstr(s);
end
