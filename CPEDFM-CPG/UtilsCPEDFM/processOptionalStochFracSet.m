function [accepted, nonPlanar, fracArea] = processOptionalStochFracSet(fracSet)
%PROCESSOPTIONALSTOCHFRACSET Safely process an optional ADFNE fracture set.
%
% Eagle Ford examples use three independently generated stochastic sets.
% ADFNE can reject/abort one set after earlier sets have succeeded. Missing
% sets represent zero optional fractures and must not abort the simulation.

if isempty(fracSet)
    accepted = [];
    nonPlanar = [];
    fracArea = zeros(0, 1);
    return
end

[accepted, nonPlanar, fracArea] = processStochFracs(fracSet);
end
