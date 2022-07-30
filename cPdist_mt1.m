clear;close all;

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

base_path = 'G:/PNAS/';

Gs = cell(2,1);

dire = dir(fullfile(base_path, '*.mat'));


for j = 1:2
    Gs{j} = load(fullfile(base_path,[dire(j).name]));
    Gs{j} = Gs{j}.G;
end
tic;rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},options);toc;
tic;rslt21 = Gs{2}.ComputeContinuousProcrustes(Gs{1},options);toc;

