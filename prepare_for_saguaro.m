clear
close all;

base_path = 'G:\ABIDE_I';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\mt1\meshes';
%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\radius\meshes';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC14\Real\Data';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\GPLmkBDMatch-master\sample_data';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';
start_path = 'G:\ABIDE_I';
% start_path = 'H:\ADNI2_01_20';
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
%dire = dir(base_path);
%dire(1:2) = [];
numLmk = 1500;
% Lmk_needed = 60;
BNN = 200;
params.preprocess.dataFileType = 'off';
% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

for i = 1:numberOfFolders
    %outputfilename = [fullfile(base_path,[dire(i).name]),'_Fan_landmark.txt'];
    dir = SubFolderName{i};
    new_dir = fullfile('G:\ABIDE_I\ABIDE_saguaro\ASD',dire(i).name);
    mkdir(new_dir);
    source = fullfile(dir,'lh_combine_final2.1.node');
    target = new_dir;
    copyfile(source, target);
    source = fullfile(dir,'lh_combine_final2.1.ele');
    copyfile(source, target);
    source = fullfile(dir,'lh_combine_final2.1.face');
    copyfile(source, target);
end



