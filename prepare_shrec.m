clear;close all;

base_path = 'G:\SHREC14\Real\Data\';
dire = dir([base_path, '*.off']);

% SIHKS
lmk = [20, 50, 80, 150, 250, 350, 450];
method = 'fan';%'fan';%, 'gao', 'euclidean', 'sm', 'ms', 'nc'};
combos = combntns(1:10,2);
num_points = 450;
%category = load(fullfile(base_path,[method,'_category_',int2str(num_points),'.mat']));
category = load(fullfile(base_path,['category_',int2str(num_points),'.mat']));
%sihks = category.category{1,1}.sihks;
%wks = category.category{1,1}.wks;
acc = 0;
sen = 0;
spe = 0;
GP_acc = 0;
feature=[];
species = [];
for k = 1:10
    class = category.category{k}.wks;
%     class = class;
%     [pc,score,latent,tsquare] = pca(class);
%     est = cumsum(latent)./sum(latent);
%     tran=pc(:,1:39);
%     class= bsxfun(@minus,class,mean(class,1));
%     class= class*tran;
    feature = [feature; class];
    species = [species; ones(40,1) * k];
end
%for l = 1:length(lmk)
    for k = 1:size(combos,1)
        class1 = category.fan_category{combos(k, 1)}.wks;
        class2 = category.fan_category{combos(k, 2)}.wks;
        feature = [feature; class1 ; class2];
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
            [out1, out2, out3, out4] = binaryEPGP([3.0; 0.0], 'covSEiso', feature(idxTrain,:), species(idxTrain), feature(idxTest,:));
            GP_acc = GP_acc + out1;
        end
        acc = acc + cp.CorrectRate;
        sen = sen + cp.Sensitivity;
        spe = spe + cp.Specificity;
    end
%end
acc = acc / 45;
sen = sen / 45;
spe = spe / 45;
GP_acc = GP_acc/450;

