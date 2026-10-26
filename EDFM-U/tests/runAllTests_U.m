%% runAllTests_U
% Run the whole EDFM-U test suite and report a single verdict.
%
%   run modules/EDFM-U/tests/runAllTests_U.m
%
% Each test is a standalone script that raises an error on failure; this
% harness catches that so one failure does not hide the rest.
%
% The tests open with `clear`, and a script invoked by `run` executes in
% its caller's workspace - so calling them directly from here would wipe
% this harness's own variables mid-loop. Each one is therefore launched
% inside runOne(), whose workspace is the one that gets cleared.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));

tests = {'test01_geometry_U', ...
         'test02_cartesian_crosscheck_U', ...
         'test03_davg_accuracy_U', ...
         'test04_pebi_dfm_validation_U'};

results = cell(numel(tests), 1);
timings = zeros(numel(tests), 1);

for i = 1:numel(tests)
    fprintf('\n########## %s ##########\n', tests{i});
    t0 = tic;
    try
        runOne(fullfile(thisDir, [tests{i}, '.m']));
        results{i} = 'PASS';
    catch me
        results{i} = ['FAIL: ', me.message];
    end
    timings(i) = toc(t0);
    close all
end

fprintf('\n================= EDFM-U TEST SUITE =================\n');
nFail = 0;
for i = 1:numel(tests)
    fprintf('  %-34s %6.1f s   %s\n', tests{i}, timings(i), results{i});
    nFail = nFail + ~strcmp(results{i}, 'PASS');
end
if nFail == 0
    fprintf('=== ALL %d TESTS PASSED (%.1f s total) ===\n', numel(tests), sum(timings));
else
    fprintf('=== %d of %d TESTS FAILED ===\n', nFail, numel(tests));
    error('EDFM-U test suite failed');
end

%% ------------------------------------------------------------------------
function runOne(scriptPath)
% Sacrificial workspace: the test's `clear` empties this one, not the
% harness's. Nothing is needed here after the script returns.
run(scriptPath);
end
