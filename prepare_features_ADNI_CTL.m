clear;close all;

AD_surface_WKS = [];
AD_tet_WKS = [];
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

for n = 1:numberOfFolders
    if exist([SubFolderName{n},'\WKS_white.txt'],'file')&&exist([SubFolderName{n},'\WKS.txt'],'file')
        dir = SubFolderName{n};
        display(SubFolderName{n});
        %cd(dir);
        pial_WKS = load([dir,'\WKS_pial.txt']);
        white_WKS = load([dir,'\WKS_white.txt']);
        volume_WKS = load([dir,'\WKS.txt']);
        
        surface_WKS = [reshape(pca(pial_WKS),1,[]), reshape(pca(white_WKS),1,[])];
        tet_WKS = reshape(pca(volume_WKS),1,[]);
        AD_surface_WKS = [AD_surface_WKS; surface_WKS];
        AD_tet_WKS = [AD_tet_WKS;tet_WKS];
    end
end
parsave([start_path,'\CTL_surface_WKS.txt'],AD_surface_WKS);
parsave([start_path,'\CTL_tet_WKS.txt'],AD_tet_WKS);
