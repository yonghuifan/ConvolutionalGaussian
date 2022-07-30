clear;close all;

base_path = 'G:\ADNI2_processed';
%dire = dir([base_path, '*.off']);

% SIHKS

num_points = 1000;
category1 = load(fullfile(base_path,'AD_.txt'));
category2 = load(fullfile(base_path,'CTL_WKS.txt'));

acc = 0;
sen = 0;
spe = 0;
%for l = 1:length(lmk)
%for k = 1:size(combos,1)
class1 = category1';
class2 = category2';

feature = [class1 ; class2];
[sample_size, feature_dim] = size(feature);
species = zeros(sample_size,1);
for i = 1:sample_size
    if i <= size(class1,1)
        species(i) = 1;
    else
        species(i) = 0;
    end
end
k=10;
cvFolds = crossvalind('Kfold', species, k);   %# get indices of 10-fold CV
cp = classperf(species);                      %# init performance tracker
rng default
cvp = cvpartition(size(feature,1),'kfold',k);
L = 0;
corL = 0;
for i = 1:k
    Train = training(cvp,i); % Extract training set indices
    [idxTrain,col] = find(Train);
    
    %     SVMModel = fitclinear(feature(idxTrain,:), species(idxTrain),'Solver','sparsa',...
    %         'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    %         struct('AcquisitionFunctionName','expected-improvement-plus'));
    SVMModel = fitcsvm(feature(idxTrain,:),species(idxTrain));
    %             SVMModel = fitcsvm(feature(idxTrain,:),species(idxTrain),'OptimizeHyperparameters','auto',...
    %                 'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    %                 'expected-improvement-plus'));
    
    CVSVMModel = crossval(SVMModel);
    classLoss = 1-kfoldLoss(CVSVMModel);
    Test = test(cvp,i);
    [idxTest,col] = find(Test);
    labels = predict(SVMModel,feature(idxTest,:));
    L = loss(SVMModel,feature(idxTest,:),species(idxTest));
    classperf(cp,labels,idxTest);
    acc = acc + cp.CorrectRate;
    sen = sen + cp.Sensitivity;
    spe = spe + cp.Specificity;
end

%    end
%end
acc = acc / 10;
sen = sen / 10;
spe = spe / 10;


