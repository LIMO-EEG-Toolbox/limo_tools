function tests = test_command_line
% Run with EEGLAB and its LIMO dependencies initialized on the MATLAB path.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'external')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'external', 'psom')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'limo_cluster_functions')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'tests', 'fixtures', 'serial_batch')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'tests', 'fixtures', 'no_dialogs')));
% Use a session override, never change the user's parallel preferences.
matlabSettings = settings;
poolSetting = matlabSettings.parallel.client.pool.AutoCreate;
hadTemporary = hasTemporaryValue(poolSetting);
oldTemporary = [];
if hadTemporary, oldTemporary = poolSetting.TemporaryValue; end
testCase.addTeardown(@restorePoolSetting, poolSetting, hadTemporary, oldTemporary);
poolSetting.TemporaryValue = false;
end

function restorePoolSetting(setting, existed, value)
if existed
    setting.TemporaryValue = value;
else
    clearTemporaryValue(setting);
end
end

function setup(testCase)
testCase.TestData.folder = pwd;
testCase.TestData.random = rng;
testCase.TestData.warnings = warning;
testCase.TestData.figures = findall(groot, 'Type', 'figure');
testCase.TestData.visibility = get(groot, 'DefaultFigureVisible');
testCase.TestData.hadStudy = evalin('base', 'exist(''STUDY'', ''var'') == 1');
if testCase.TestData.hadStudy
    testCase.TestData.study = evalin('base', 'STUDY');
end
evalin('base', 'clear STUDY');
testCase.TestData.scratch = tempname;
mkdir(testCase.TestData.scratch);
cd(testCase.TestData.scratch);
set(groot, 'DefaultFigureVisible', 'off');
rng(42);
chanlocs = struct('labels', {'C1', 'C2'}, 'X', {1, -1}, 'Y', {0, 0}, 'Z', {1, 1});
testCase.TestData.channels = struct('expected_chanlocs', chanlocs, 'channeighbstructmat', [0 1; 1 0]);
for subject = 1:18
    folder = fullfile(testCase.TestData.scratch, sprintf('subject_%02d', subject));
    mkdir(folder);
    LIMO = struct('Analysis', 'Time', 'Type', 'Channels', 'Level', 1, 'dir', folder);
    LIMO.data = struct('sampling_rate', 100, 'trim1', 1, 'trim2', 5, 'start', 0, 'end', 40, 'chanlocs', chanlocs);
    LIMO.design.X = [randn(20, 6) ones(20, 1)];
    LIMO.design.labels = arrayfun(@(n) struct('description', sprintf('parameter %d', n)), 1:7);
    save(fullfile(folder, 'LIMO.mat'), 'LIMO');
    Betas = randn(2, 5, 7) + subject / 10;
    save(fullfile(folder, 'Betas.mat'), 'Betas');
    testCase.TestData.betas{subject} = fullfile(folder, 'Betas.mat');
    for parameter = 1:2
        con = cat(3, Betas(:, :, parameter), zeros(2, 5, 4));
        file = fullfile(folder, sprintf('con_%d.mat', parameter));
        save(file, 'con');
        testCase.TestData.con{subject, parameter} = file;
    end
end
end

function teardown(testCase)
cd(testCase.TestData.folder);
delete(setdiff(findall(groot, 'Type', 'figure'), testCase.TestData.figures));
set(groot, 'DefaultFigureVisible', testCase.TestData.visibility);
rng(testCase.TestData.random);
warning(testCase.TestData.warnings);
evalin('base', 'clear STUDY');
if testCase.TestData.hadStudy
    assignin('base', 'STUDY', testCase.TestData.study);
end
rmdir(testCase.TestData.scratch, 's');
end

function testSettingsPreserveExplicitStudyWithEmptyBase(testCase)
assignin('base', 'STUDY', []);
STUDY = struct('filepath', pwd, 'filename', 'explicit.study');
expected = STUDY;
limo_settings_script;
verifyEqual(testCase, STUDY, expected);
verifyEqual(testCase, limo_settings.workdir, fullfile(pwd, 'derivatives'));
end

function testSettingsPreserveExplicitStudyWithDifferentBase(testCase)
assignin('base', 'STUDY', struct('filepath', tempdir, 'filename', 'other.study'));
STUDY = struct('filepath', pwd, 'filename', 'explicit.study');
expected = STUDY;
limo_settings_script;
verifyEqual(testCase, STUDY, expected);
verifyEqual(testCase, limo_settings.workdir, fullfile(pwd, 'derivatives'));
end

function testSettingsDoNotCreateEmptyStudy(testCase)
assignin('base', 'STUDY', []);
limo_settings_script;
verifyEqual(testCase, exist('STUDY', 'var'), 0);
verifyEqual(testCase, limo_settings.workdir, '');
end

function testCustomSettingsPreserveExplicitStudy(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'tests', 'fixtures', 'custom_settings')));
STUDY = struct('filepath', pwd, 'filename', 'explicit.study');
expected = STUDY;
limo_settings_script;
verifyEqual(testCase, STUDY, expected);
verifyEqual(testCase, limo_settings.workdir, fullfile(tempdir, 'limo_custom_output'));
end

function testStandaloneContrastsWithoutStudy(testCase)
runContrast(testCase);
end

function testFirstLevelSubjectFilename(testCase)
folder = fullfile(testCase.TestData.scratch, 'sub-001', 'eeg');
mkdir(folder);
EEG = eeg_emptyset;
EEG.data = randn(2, 5, 40);
EEG.nbchan = 2; EEG.pnts = 5; EEG.trials = 40;
EEG.srate = 100; EEG.xmin = 0; EEG.xmax = 0.04; EEG.times = 0:10:40;
EEG.chanlocs = testCase.TestData.channels.expected_chanlocs;
EEG.filename = 'sub-001_task-test_eeg.set';
EEG.filepath = folder;
EEG.etc.timeerp = EEG.times;
EEG.etc.datafiles.daterp = fullfile(folder, 'signal.mat');
signal = EEG.data;
save(EEG.etc.datafiles.daterp, 'signal');
file = fullfile(folder, EEG.filename);
save(file, 'EEG');
STUDY = struct('filepath', testCase.TestData.scratch, 'filename', 'test.study', ...
    'currentdesign', 1, 'design', struct('name', 'Synthetic'), 'group', {{'1'}});
STUDY.datasetinfo = struct('filename', EEG.filename, 'filepath', folder, ...
    'subject', 'sub-001', 'session', [], 'group', '1');
model = struct('set_files', {{file}}, 'cat_files', {{[ones(20,1); 2*ones(20,1)]}}, 'cont_files', []);
model.defaults = struct('analysis', 'Time', 'type', 'Channels', 'zscore', 0, ...
    'method', 'OLS', 'type_of_analysis', 'Mass-univariate', 'fullfactorial', 0, ...
    'bootstrap', 0, 'tfce', 0, 'start', 0, 'end', 40, 'verbose', 'noGUI');
assignin('base', 'STUDY', []);
[files, status] = limo_batch('model specification', model, [], STUDY);
verifyEqual(testCase, status, 1);
verifyTrue(testCase, isfile(files.Beta{1}));
[~, name, ext] = fileparts(files.Beta{1});
verifyEqual(testCase, [name ext], 'sub-001_desc-Betas.mat');
end

function testStandaloneContrastsWithEmptyStudy(testCase)
assignin('base', 'STUDY', []);
runContrast(testCase);
end

function testMissingContrastMatrixErrorsWithoutDialog(testCase)
contrast = struct('LIMO_files', {{fullfile(pwd, 'LIMO.mat')}});
verifyError(testCase, @() limo_batch('contrast only', [], contrast), 'LIMO:InvalidContrast');
end

function testBestElectrodesPrefixedResults(testCase)
checkBestElectrodes(testCase, true);
end

function testBestElectrodesLegacyResults(testCase)
checkBestElectrodes(testCase, false);
end

function checkBestElectrodes(testCase, prefixed)
folder = fullfile(testCase.TestData.scratch, 'sub-001', 'eeg', 'model');
mkdir(folder);
loaded = load(fullfile(fileparts(testCase.TestData.betas{1}), 'LIMO.mat'));
LIMO = loaded.LIMO;
LIMO.dir = folder;
file = fullfile(folder, 'LIMO.mat');
save(file, 'LIMO');
R2 = zeros(2, 5, 3);
R2(2,3,2) = 10;
if prefixed
    save(fullfile(folder, 'sub-001_desc-R2.mat'), 'R2');
    R2 = flip(R2, 1); % stale legacy output must not take precedence
end
save(fullfile(folder, 'R2.mat'), 'R2');
[channel, original] = limo_best_electrodes({file});
verifyEqual(testCase, channel, 2);
verifyTrue(testCase, isnan(original));
end

function runContrast(testCase)
folder = fullfile(testCase.TestData.scratch, 'sub-001', 'eeg', 'model');
mkdir(folder);
loaded = load(fullfile(fileparts(testCase.TestData.betas{1}), 'LIMO.mat'));
LIMO = loaded.LIMO;
LIMO.dir = folder;
LIMO.design.name = 'Synthetic OLS';
LIMO.design.method = 'OLS';
LIMO.design.type_of_analysis = 'Mass-univariate';
LIMO.design.bootstrap = 0;
LIMO.design.tfce = 0;
LIMO.model.model_df = repmat([7 13], 2, 1);
loaded = load(testCase.TestData.betas{1});
Betas = loaded.Betas;
Res = randn(2, 5, 20);
Yr = Res;
for channel = 1:2
    Yr(channel, :, :) = squeeze(Res(channel, :, :)) + squeeze(Betas(channel, :, :)) * LIMO.design.X';
end
file = fullfile(folder, 'LIMO.mat');
save(file, 'LIMO');
save(fullfile(folder, 'sub-001_desc-Yr.mat'), 'Yr');
save(fullfile(folder, 'sub-001_desc-Betas.mat'), 'Betas');
save(fullfile(folder, 'sub-001_desc-Res.mat'), 'Res');
contrast = struct('LIMO_files', {{file}}, 'mat', [1 -1 0 0 0 0 0]);
[files, status] = limo_batch('contrast only', [], contrast);
verifyEqual(testCase, status, 1);
verifyTrue(testCase, isfile(files.con{1}{1}));
end

function testOneSampleMatchesTrimmedMeanReference(testCase)
folder = runAnalysis(testCase, 'one sample t-test', testCase.TestData.con(:, 1));
result = load(fullfile(folder, 'One_Sample_Ttest_parameter_1.mat'));
for channel = 1:2
    for frame = 1:5
        values = zeros(1, 18);
        for subject = 1:18
            data = load(testCase.TestData.con{subject, 1});
            values(subject) = data.con(channel, frame, 1);
        end
        values = sort(values);
        trim = floor(0.2 * numel(values));
        retained = numel(values) - 2 * trim;
        center = mean(values(trim + 1:end - trim));
        winsorized = values;
        winsorized(1:trim) = values(trim + 1);
        winsorized(end - trim + 1:end) = values(end - trim);
        % Preserve LIMO's existing asymptotic standard error convention.
        se = std(winsorized) / (0.6 * sqrt(numel(values)));
        statistic = center / se;
        expected = [center, se, retained - 1, statistic, 2 * tcdf(-abs(statistic), retained - 1)];
        verifyEqual(testCase, reshape(result.one_sample(channel, frame, :), 1, []), expected, 'AbsTol', 1e-10);
    end
end
end

function testPairedNestedScalarPaths(testCase)
files = cellfun(@(file) {file}, testCase.TestData.con, 'UniformOutput', false);
folder = runAnalysis(testCase, 'paired t-test', files);
verifyNotEmpty(testCase, dir(fullfile(folder, 'Paired_Samples_Ttest*.mat')));
end

function testPairedMissingPartnersAreNotRealigned(testCase)
files = testCase.TestData.con;
files{2,1} = '';
files{3,2} = {''};
verifyError(testCase, @() runAnalysis(testCase, 'paired t-test', files), 'LIMO:UnmatchedPairs');
end

function testPairedEmptyRowsPreserveAlignment(testCase)
files = testCase.TestData.con;
files(2,:) = {''};
folder = runAnalysis(testCase, 'paired t-test', files);
first = load(fullfile(folder, 'Y1r.mat'));
second = load(fullfile(folder, 'Y2r.mat'));
expected = zeros(2, 5, 17);
subjects = [1 3:18];
for index = 1:numel(subjects)
    a = load(files{subjects(index),1});
    b = load(files{subjects(index),2});
    expected(:,:,index) = a.con(:,:,1) - b.con(:,:,1);
end
verifyEqual(testCase, first.Y1r - second.Y2r, expected, 'AbsTol', 1e-12);
end

function testTwoSamples(testCase)
files = reshape(testCase.TestData.con(:, 1), 9, 2)';
folder = runAnalysis(testCase, 'two-samples t-test', files);
verifyNotEmpty(testCase, dir(fullfile(folder, 'Two_Samples_Ttest*.mat')));
end

function testRegression(testCase)
folder = runAnalysis(testCase, 'regression', testCase.TestData.betas, ...
    'parameters', 1, 'regressor', (1:18)');
verifyTrue(testCase, isfile(fullfile(folder, 'R2.mat')));
actual = load(fullfile(folder, 'Yr.mat'));
expected = zeros(2, 5, 18);
for subject = 1:18
    loaded = load(testCase.TestData.betas{subject});
    expected(:,:,subject) = loaded.Betas(:,:,1);
end
verifyEqual(testCase, actual.Yr, expected, 'AbsTol', 1e-12);
end

function testOptimizedNumericChannels(testCase)
channels = repmat([1; 2], 9, 1);
folder = runAnalysis(testCase, 'regression', testCase.TestData.betas, ...
    'parameters', 1, 'regressor', (1:18)', ...
    'analysis_type', '1 channel/component only', 'channel', channels);
verifySelectedChannels(testCase, folder, channels, testCase.TestData.betas);
end

function testOptimizedChannelFile(testCase)
channels = repmat([1 2], 1, 9);
file = fullfile(pwd, 'channels.mat');
save(file, 'channels');
folder = runAnalysis(testCase, 'one sample t-test', testCase.TestData.con(:,1), ...
    'analysis_type', '1 channel/component only', 'channel', file);
verifySelectedChannels(testCase, folder, channels, testCase.TestData.con(:,1));
end

function testOptimizedGroupedChannelsPreserveOrder(testCase)
channels = [ones(1,7) 2*ones(1,11)];
files = cell(2,11);
files(1,1:7) = testCase.TestData.con(1:7,1)';
files(2,:) = testCase.TestData.con(8:18,1)';
folder = runAnalysis(testCase, 'two-samples t-test', files, ...
    'analysis_type', '1 channel/component only', 'channel', channels);
first = load(fullfile(folder, 'Y1r.mat'));
second = load(fullfile(folder, 'Y2r.mat'));
expected = zeros(1,5,18);
for subject = 1:18
    loaded = load(testCase.TestData.con{subject,1});
    expected(1,:,subject) = loaded.con(channels(subject),:,1);
end
verifyEqual(testCase, cat(3,first.Y1r,second.Y2r), expected, 'AbsTol', 1e-12);
end

function verifySelectedChannels(testCase, folder, channels, files)
actual = load(fullfile(folder, 'Yr.mat'));
expected = zeros(1, 5, 18);
for subject = 1:18
    loaded = load(files{subject});
    values = loaded.(cell2mat(fieldnames(loaded)));
    expected(1,:,subject) = values(channels(subject),:,1);
end
verifyEqual(testCase, actual.Yr, expected, 'AbsTol', 1e-12);
end

function testSingleChannelTimeTfceBootstrap(testCase)
verifySingleChannelTfce(testCase, false);
end

function testCentralSummaryReturnsWithoutPopups(testCase)
blockSummaryGui(testCase);
data = randn(2,5,18);
result = limo_central_tendency_and_ci(data, 'Mean', []);
verifyEqual(testCase, result.mean(:,:,1,2), mean(data,3), 'AbsTol', 1e-12);
verifyEqual(testCase, findall(groot, 'Type', 'figure'), testCase.TestData.figures);
end

function testCentralSummarySavesWithoutPopups(testCase)
blockSummaryGui(testCase);
data = randn(2,5,18);
limo_central_tendency_and_ci(data, 'Mean', [], 'summary');
loaded = load('summary_Mean.mat');
verifyEqual(testCase, loaded.Data.mean(:,:,1,2), mean(data,3), 'AbsTol', 1e-12);
verifyEqual(testCase, findall(groot, 'Type', 'figure'), testCase.TestData.figures);
end

function testCentralSummaryTimeFrequencyWithoutPopups(testCase)
blockSummaryGui(testCase);
data = randn(2,3,5,18);
result = limo_central_tendency_and_ci(data, 'Mean', []);
verifyEqual(testCase, result.mean(:,:,:,1,2), mean(data,4), 'AbsTol', 1e-12);
verifyEqual(testCase, findall(groot, 'Type', 'figure'), testCase.TestData.figures);
end

function blockSummaryGui(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'tests','fixtures','no_summary_gui')));
end

function testPlotSuppliedFilesReturnsWithoutChooser(testCase)
checkPlotSelection(testCase, 'mean', 1);
end

function testPlotSuppliedVariableDoesNotReplaceFiles(testCase)
checkPlotSelection(testCase, 'mean', 2);
end

function testPlotSuppliedSubjectSelectsSubjectAxis(testCase)
checkPlotSelection(testCase, 'data', 3);
end

function checkPlotSelection(testCase, field, variable)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'tests','fixtures','no_summary_gui')));
LIMO = struct('Analysis','Time','data',struct('start',0,'end',40,'sampling_rate',100));
files = cell(1,2);
expected = zeros(2,5);
for index = 1:2
    Data.limo = LIMO;
    if strcmp(field,'mean')
        Data.mean = zeros(2,5,4,3);
        for condition = 1:4
            center = (1:5) + index + condition*10;
            Data.mean(2,:,condition,1) = center - 1;
            Data.mean(2,:,condition,2) = center;
            Data.mean(2,:,condition,3) = center + 1;
        end
        expected(index,:) = Data.mean(2,:,variable,2);
    else
        Data.data = reshape(1:40,2,5,1,4) + index;
        expected(index,:) = Data.data(2,:,1,variable);
    end
    files{index} = fullfile(pwd,sprintf('summary_%d.mat',index));
    save(files{index},'Data');
end
limo_add_plots(files,LIMO,'channel',2,'variable',variable);
figures = setdiff(findall(groot,'Type','figure'),testCase.TestData.figures);
verifyNumElements(testCase,figures,1);
lines = findall(figures,'Type','line');
verifyNumElements(testCase,lines,2);
actual = get(lines,'YData');
verifyEqual(testCase,sortrows(vertcat(actual{:})),sortrows(expected),'AbsTol',1e-12);
end

function testSingleChannelTimeFrequencyTfceBootstrap(testCase)
verifySingleChannelTfce(testCase, true);
end

function verifySingleChannelTfce(testCase, timeFrequency)
LIMO = struct('Level', 2, 'dir', pwd, 'Analysis', 'Time');
LIMO.design.name = 'One sample t-test';
LIMO.data.neighbouring_matrix = 0;
if timeFrequency
    LIMO.Analysis = 'Time-Frequency';
    one_sample = zeros(1, 3, 5, 5);
    one_sample(:,:,:,4) = reshape(1:15, 1, 3, 5) / 5;
    H0_one_sample = zeros(1, 3, 5, 2, 3);
    for bootstrap = 1:3
        H0_one_sample(:,:,:,1,bootstrap) = one_sample(:,:,:,4) / bootstrap;
    end
else
    one_sample = zeros(1, 5, 5);
    one_sample(:,:,4) = [1 2 3 2 1];
    H0_one_sample = zeros(1, 5, 2, 3);
    for bootstrap = 1:3
        H0_one_sample(:,:,1,bootstrap) = one_sample(:,:,4) / bootstrap;
    end
end
save('LIMO.mat', 'LIMO');
mkdir('H0');
name = 'One_Sample_Ttest_parameter_1';
save([name '.mat'], 'one_sample');
save(fullfile('H0', [name '_desc-H0.mat']), 'H0_one_sample');
limo_tfce_handling([name '.mat']);
actual = load(fullfile('H0', [name '_desc-tfceH0.mat']));
if timeFrequency
    verifySize(testCase, actual.tfce_H0_score, [1 3 5 3]);
else
    verifySize(testCase, actual.tfce_H0_score, [1 5 3]);
end
for bootstrap = 1:3
    if timeFrequency
        expected = limo_tfce(2, squeeze(H0_one_sample(:,:,:,1,bootstrap)), [], 0);
        verifyEqual(testCase, squeeze(actual.tfce_H0_score(:,:,:,bootstrap)), expected, 'AbsTol', 1e-12);
    else
        expected = limo_tfce(1, H0_one_sample(:,:,1,bootstrap), 0, 0);
        verifyEqual(testCase, actual.tfce_H0_score(:,:,bootstrap), expected, 'AbsTol', 1e-12);
    end
end
end

function testAnovaFileLists(testCase)
files = cell(1, 3);
for group = 1:3
    files{group} = fullfile(pwd, sprintf('group_%d.txt', group));
    writelines(string(testCase.TestData.con((group - 1) * 6 + (1:6), 1)), files{group});
end
folder = runAnalysis(testCase, 'N-Ways ANOVA', files);
verifyTrue(testCase, isfile(fullfile(folder, 'LIMO.mat')));
end

function testAncovaRaggedNestedPaths(testCase)
files = cell(3, 7);
groups = {1:6, 7:13, 14:18};
for group = 1:3
    for subject = 1:numel(groups{group})
        files{group, subject} = testCase.TestData.con(groups{group}(subject), 1);
    end
end
folder = runAnalysis(testCase, 'ANCOVA', files, 'regressor', (1:18)');
verifyTrue(testCase, isfile(fullfile(folder, 'R2.mat')));
end

function testRepeatedMeasuresUsesSuppliedParameters(testCase)
folder = runAnalysis(testCase, 'Repeated Measures ANOVA', testCase.TestData.betas, ...
    'parameters', {[1 2 3], [4 5 6]}, 'factor names', {'condition', 'repetition'});
verifyNotEmpty(testCase, dir(fullfile(folder, 'Rep_ANOVA*.mat')));
end

function testRepeatedMeasuresFileList(testCase)
file = fullfile(pwd, 'Beta_files.txt');
writelines(string(testCase.TestData.betas), file);
folder = runAnalysis(testCase, 'Repeated Measures ANOVA', {file}, ...
    'parameters', {[1 2 3], [4 5 6]}, 'factor names', {'condition', 'repetition'});
verifyNotEmpty(testCase, dir(fullfile(folder, 'Rep_ANOVA*.mat')));
end

function testRepeatedMeasuresGroupedFileLists(testCase)
files = cell(3, 1);
for group = 1:3
    files{group} = fullfile(pwd, sprintf('Beta_group_%d.txt', group));
    writelines(string(testCase.TestData.betas((group-1)*6+(1:6))), files{group});
end
folder = runAnalysis(testCase, 'Repeated Measures ANOVA', files, ...
    'parameters', {[1 2 3]}, 'factor names', {'condition'});
verifyTrue(testCase, isfile(fullfile(folder, 'Rep_ANOVA_Gp_effect.mat')));
verifyRepeatedMeasuresOrder(testCase, folder, 1:2);
end

function testRepeatedMeasuresOptimizedGroupedChannels(testCase)
files = cell(3,1);
groups = {1:6, 7:13, 14:18};
for group = 1:3
    files{group} = fullfile(pwd, sprintf('Beta_group_%d.txt', group));
    writelines(string(testCase.TestData.betas(groups{group})), files{group});
end
channels = [ones(1,6) 2*ones(1,7) ones(1,5)];
folder = runAnalysis(testCase, 'Repeated Measures ANOVA', files, ...
    'parameters', {[1 2 3]}, 'factor names', {'condition'}, ...
    'analysis_type', '1 channel/component only', 'channel', channels);
verifyRepeatedMeasuresOrder(testCase, folder, channels);
end

function verifyRepeatedMeasuresOrder(testCase, folder, channels)
actual = load(fullfile(folder, 'Yr.mat'));
if numel(channels) == 18
    expected = zeros(1,5,18,3);
else
    expected = zeros(2,5,18,3);
end
for subject = 1:18
    loaded = load(testCase.TestData.betas{subject});
    if numel(channels) == 18
        selected = channels(subject);
    else
        selected = channels;
    end
    expected(:,:,subject,:) = loaded.Betas(selected,:,1:3);
end
verifyEqual(testCase, actual.Yr, expected, 'AbsTol', 1e-12);
end

function folder = runAnalysis(testCase, kind, files, varargin)
folder = fullfile(testCase.TestData.scratch, 'analysis');
mkdir(folder);
cd(folder);
result = limo_random_select(kind, testCase.TestData.channels, 'LIMOfiles', files, ...
    'analysis_type', 'Full scalp analysis', 'type', 'Channels', 'nboot', 0, ...
    'tfce', 0, 'skip design check', 'yes', 'zscore', 'no', varargin{:});
verifyNotEmpty(testCase, result);
end
