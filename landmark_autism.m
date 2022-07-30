%%% init
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
    if ~exist(fullfile(dir,'_fan_landmark2000.txt'),'file')
        display(dire(i).name);
        filename = fullfile(dir,'lh_combine_final2.1');
        [node,elem,~]=readtetgen(filename);
        face = [];
%         parfor j = 1:size(elem,1)
%             face = [face;sort([elem(j,1),elem(j, 2), elem(j,3)])];
%             face = [face;sort([elem(j,1),elem(j, 2), elem(j,4)])];
%             face = [face;sort([elem(j,1),elem(j, 3), elem(j,4)])];
%             face = [face;sort([elem(j,4),elem(j, 2), elem(j,3)])];
%         end
%         face = unique(face,'rows');
        saveoff(node,face,fullfile(dir,'node.off'));
    %end
         G = Mesh(params.preprocess.dataFileType,fullfile(dir,'node.off'));
         node = G.V';
    %
         [idx, dist] = knnsearch(node,node,'K',BNN+1);
    %
         Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
    %     % Lmks = GetGPLmk(G, numLmk);
    %     % Lmks = GetGPLmk_Euclidean(G, numLmk);
    %     %Lmks = GetGPLmk_NoCurvature(G, numLmk);
    %     % Lmks = GetGPLmk_NoCurvature(G, numLmk);
    %     Lmks = reshape(Lmks,1,[]);
    %
    save(fullfile(dir,'_fan_landmark2000.txt'),'Lmks','-ascii');
end
end
