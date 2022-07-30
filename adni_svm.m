clear;close all;

base_path = pwd;
dire = dir([base_path, '*.off']);

% SIHKS
lmk = [20, 50, 80, 150, 250, 350, 450];
method = '';%'fan';%, 'gao', 'euclidean', 'sm', 'ms', 'nc'};

num_points = 1000;
category1 = load(fullfile(pwd,['AD_',int2str(num_points),'.mat']));
category2 = load(fullfile(pwd,['CTL_',int2str(num_points),'.mat']));
%category = load(fullfile(base_path,[method,'category_',int2str(num_points),'.mat']));
%sihks = category.category{1,1}.sihks;
%wks = category.category{1,1}.wks;
acc = 0;
sen = 0;
spe = 0;
%for l = 1:length(lmk)
%for k = 1:size(combos,1)
class1 = category1.wks_tmp;
class2 = category2.wks_tmp;
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


