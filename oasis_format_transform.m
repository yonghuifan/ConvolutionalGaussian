clear;close all;

% start_path = 'G:\Amyloid_matched\NL_pos';
start_path = 'G:\OASIS\pos';
toplevel = uigetdir(start_path);

dire = dir(toplevel);
%dire(~[dire.isdir]) = [];
dire = dire(3:end);

% subFolder = genpath(toplevel);
% remain = subFolder;
% listOfFolderNames = {};
% 
% while(true)
%     [singleSubFolder, remain] = strtok(remain, ';');
%     if isempty(singleSubFolder)
%         break;
%     end
%     listOfFolderNames = [listOfFolderNames singleSubFolder];
% end
% j = 1;
% for i = 2:1:length(listOfFolderNames)
%     SubFolderName{j} = listOfFolderNames{i};
%     j = j+1;
% end
% SubFolderName = dire(:).name;
numberOfFolders = length(dire);
% Guadalupe_dir = [userpath,'\Guadalupe2.exe'];
for n = 1:numberOfFolders
    tic
    dir = dire(n).name;
    display(dir);
    SubFolderName = [toplevel,'\',dir];
    output_name = [SubFolderName,'.off'];
    oasis_m_to_off(SubFolderName, output_name);
    toc
end