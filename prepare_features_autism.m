clear;close all;

ASD_WKS = [];
TC_WKS = [];
start_path = 'G:\ABIDE_I';
ASD_path1 = 'G:\ABIDE_I\ASD\';
ASD_path2 = 'G:\ABIDE_I\ASD_2\';
TC_path1 = 'G:\ABIDE_I\TC\';
TC_path2 = 'G:\ABIDE_I\TC_2\';
dire = dir(ASD_path1);
dire(~[dire.isdir]) = [];
dire = dire(3:end);

params.preprocess.dataFileType = 'off';

parfor n = 1:length(dire)
    path = fullfile(ASD_path1, dire(n).name);
    if exist(fullfile(path,'\WKS_50_100.txt'),'file')&&exist(fullfile(path,'fan_landmark1500.txt'),'file')
        display(dire(n).name);
        %cd(dir);
        WKS = load(fullfile(path,'WKS_50_100.txt'));
        landmark = load(fullfile(path,'fan_landmark1500.txt'));
        WKS = WKS(landmark(1,1:1000),:);
        WKS = reshape(WKS,1,[]);
        ASD_WKS = [ASD_WKS;WKS];
    end
end
dire = dir(ASD_path2);
dire(~[dire.isdir]) = [];
dire = dire(3:end);
parfor n = 1:length(dire)
    path = fullfile(ASD_path2, dire(n).name);
    if exist(fullfile(path,'\WKS_50_100.txt'),'file')&&exist(fullfile(path,'fan_landmark1500.txt'),'file')
        display(dire(n).name);
        %cd(dir);
        WKS = load(fullfile(path,'WKS_50_100.txt'));
        landmark = load(fullfile(path,'fan_landmark1500.txt'));
        WKS = WKS(landmark(1,1:1000),:);
        WKS = reshape(WKS,1,[]);
        ASD_WKS = [ASD_WKS;WKS];
    end
end
dire = dir(TC_path1);
dire(~[dire.isdir]) = [];
dire = dire(3:end);
parfor n = 1:length(dire)
    path = fullfile(TC_path1, dire(n).name);
    if exist(fullfile(path,'\WKS_50_100.txt'),'file')&&exist(fullfile(path,'fan_landmark1500.txt'),'file')
        display(dire(n).name);
        %cd(dir);
        WKS = load(fullfile(path,'WKS_50_100.txt'));
        landmark = load(fullfile(path,'fan_landmark1500.txt'));
        WKS = WKS(landmark(1,1:1000),:);
        WKS = reshape(WKS,1,[]);
        TC_WKS = [TC_WKS;WKS];
    end
end
dire = dir(TC_path2);
dire(~[dire.isdir]) = [];
dire = dire(3:end);
parfor n = 1:length(dire)
    path = fullfile(TC_path2, dire(n).name);
    if exist(fullfile(path,'\WKS_50_100.txt'),'file')&&exist(fullfile(path,'fan_landmark1500.txt'),'file')
        display(dire(n).name);
        %cd(dir);
        WKS = load(fullfile(path,'WKS_50_100.txt'));
        landmark = load(fullfile(path,'fan_landmark1500.txt'));
        WKS = WKS(landmark(1,1:1000),:);
        WKS = reshape(WKS,1,[]);
        TC_WKS = [TC_WKS;WKS];
    end
end
parsave([start_path,'\ASD_WKS_1000.txt'],ASD_WKS);
parsave([start_path,'\TC_WKS_1000.txt'],TC_WKS);
