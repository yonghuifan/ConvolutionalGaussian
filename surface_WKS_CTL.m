clear;close all;


start_path = 'G:\ADNI2_processed\';
toplevel = uigetdir(start_path);

dire = dir(toplevel);
dire(~[dire.isdir]) = [];
dire = dire(3:end);

subFolder = genpath(toplevel);
remain = subFolder;
listOfFolderNames = {};

while(true)
    [singleSubFolder, remain] = strtok(remain, ';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder];
end
j = 1;
for i = 2:1:length(listOfFolderNames)
    SubFolderName{j} = listOfFolderNames{i};
    j = j+1;
end

numberOfFolders = length(SubFolderName);
Guadalupe_dir = [userpath,'\Guadalupe.exe'];
params.preprocess.dataFileType = 'off';

parfor n = 1:numberOfFolders
    if exist([SubFolderName{n},'\WKS_white.txt'],'file')
        dir = SubFolderName{n};
        display(SubFolderName{n});
        %cd(dir);
        [pv, pf] = readoff([dir,'\lh_pial_fixed.off']);
        [wv, wf] = readoff([dir,'\lh_white_fixed.off']);
        saveoff(pv,pf,[dir,'\lh_pial_fixed.off']);
        saveoff(wv,wf,[dir,'\lh_white_fixed.off']);
        
        G_pial = Mesh(params.preprocess.dataFileType,[dir,'\lh_pial_fixed.off']);
        G_white = Mesh(params.preprocess.dataFileType,[dir,'\lh_white_fixed.off']);
        
        % normalize (center and scale to unit area)
        G_pial.Normalize();
        G_white.Normalize();
        % reorient (outward facing normals)
        [~,~,flip] = G_pial.ComputeNormal();
        if flip
            G_pial.F = G_pial.F([1 3 2],:);
        end
        [~,~,flip] = G_white.ComputeNormal();
        if flip
            G_white.F = G_white.F([1 3 2],:);
        end
        % compute WKS
        WKS_pial = G_pial.ComputeWKS([]);
        WKS_pial = WKS_pial./max(max(WKS_pial));
        WKS_white = G_white.ComputeWKS([]);
        WKS_white = WKS_white./max(max(WKS_white));
        
        parsave([dir,'\WKS_pial.txt'],WKS_pial);
        parsave([dir,'\WKS_white.txt'],WKS_white);
    end
end