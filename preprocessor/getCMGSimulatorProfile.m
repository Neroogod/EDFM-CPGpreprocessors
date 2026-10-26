function profile = getCMGSimulatorProfile(simulator, version)
% Return CMG-family metadata used by the embedded-fracture exporter.

if nargin < 1 || isempty(simulator), simulator = 'IMEX'; end
if nargin < 2 || isempty(version), version = 'unspecified'; end
simulator = upper(validatestring(simulator, {'IMEX', 'GEM', 'STARS'}));

profile = struct;
profile.simulator = simulator;
profile.version = char(string(version));
profile.unitSystem = 'SI';
profile.transmissibilityUnit = 'md.m';
profile.includeExtension = lower([simulator, '.inc']);

switch simulator
    case 'IMEX'
        profile.physics = 'isothermal black-oil and conventional multiphase flow';
        profile.hostModel = 'IMEX black-oil dataset';
        profile.propertyPlacement = ...
            'Place grid, null, property, and multiplier records in the IMEX grid-property section.';
        profile.connectionPlacement = ...
            'Import the NNC CSV with the IMEX explicit non-neighbor-connection facility or Builder.';
        profile.phaseNote = ...
            'Use for oil-water, oil-gas, or black-oil EDFM-family cases.';
    case 'GEM'
        profile.physics = 'isothermal compositional and EOS flow';
        profile.hostModel = 'GEM compositional dataset';
        profile.propertyPlacement = ...
            'Place grid, null, property, and multiplier records in the GEM reservoir-grid section.';
        profile.connectionPlacement = ...
            'Import the NNC CSV with the GEM explicit non-neighbor-connection facility or Builder.';
        profile.phaseNote = ...
            'Define components, EOS, and compositions in the host GEM model; this exporter supplies fracture geometry and connectivity.';
    case 'STARS'
        profile.physics = 'thermal, chemical, and reactive multiphase flow';
        profile.hostModel = 'STARS thermal dataset';
        profile.propertyPlacement = ...
            'Place grid, null, property, and multiplier records in the STARS reservoir-grid section.';
        profile.connectionPlacement = ...
            'Import the NNC CSV with the STARS explicit non-neighbor-connection facility or Builder.';
        profile.phaseNote = ...
            'Define thermal properties, components, reactions, and initial temperature in the host STARS model.';
end

profile.caveat = ['CMG connection-card syntax varies by release. The ', ...
    'simulator-tagged CSV tables are the authoritative numeric output; ', ...
    'verify or import them with Builder or the keyword manual for the selected version.'];
end
