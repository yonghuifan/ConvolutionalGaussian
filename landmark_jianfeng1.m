%%% init
clear
close all;

base_path = 'G:\AD_MCI_NC\AD_pos';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\mt1\meshes';
%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\radius\meshes';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC14\Real\Data';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\GPLmkBDMatch-master\sample_data';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

dire = dir(base_path);
dire(1:2) = [];
numLmk = 1500;
% Lmk_needed = 60;
BNN = 1000;
params.preprocess.dataFileType = 'off';
% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

for i = 1:length(dire)
%outputfilename = [fullfile(base_path,[dire(i).name]),'_SM_landmark1000.txt'];
%if ~exist(outputfilename,'file')
display(dire(i).name);
filename = [fullfile(base_path,dire(i).name),'\lh_combine_final2.1'];
outname = [fullfile(base_path,dire(i).name),'\tet2surface.off'];
[node,elem,face]=readtetgen(filename);
%     DT = triangulation(elem(), node);
%     face = [];
%     parfor j = 1:size(elem,1)
%         face = [face;sort([elem(j,1),elem(j, 2), elem(j,3)])];
%         face = [face;sort([elem(j,1),elem(j, 2), elem(j,4)])];
%         face = [face;sort([elem(j,1),elem(j, 3), elem(j,4)])];
%         face = [face;sort([elem(j,4),elem(j, 2), elem(j,3)])];
%     end
%     face = unique(face,'rows');
%     writeOFF(outname,node,face);
% %end
% end
% G = Mesh(params.preprocess.dataFileType,'lh_pial_noCC.off');
G.V = node;
G.nV = size(node,1);
%node = G.V';
%face = G.F';
[idx, dist] = knnsearch(node,node,'K',BNN+1);

Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
% Lmks = GetGPLmk(G, numLmk);
% Lmks = GetGPLmk_Euclidean(G, numLmk);
%Lmks = GetGPLmk_NoCurvature(G, numLmk);
% Lmks = GetGPLmk_NoCurvature(G, numLmk);
Lmks = reshape(Lmks,1,[]);
shape = ShapeStruct(node(Lmks,:), []);
color = zeros(numLmk,3);
% color(Lmks,1) = 1;
[ fig, S ] = RenderPointCloudFunction(shape, color);
%     set(gca, 'Layer', 'top')
%     pp= patch('Faces',face,'Vertices',node,'FaceColor','interp','EdgeColor','none');
%     %pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
%     %view(90,180);
%     view(180, -90);
%     %camroll(90)
%     ax = gca;
%     %ax.CameraPosition = [0 10 0];
%     caxis([0,400]);
%     axis off;
% %     %colormap jet
% %     %colorbar
%     hold on
%     sphereSize = 0.01;
%     c = [1 0 0];
%     drawSpheres(node(Lmks,:),sphereSize,c);
% [pn,connnum,count]=meshconn(face,size(node,1));
% 
% pial_conflict = Lmks';
% %white_conflict = pos;
% count = 1;
% while count < 10
%     for i = 1:size(pial_conflict,1)
%         pial_conflict = [pial_conflict;unique(pn{pial_conflict(i)}')];
%     end
%     pial_conflict = unique(pial_conflict);
%     count = count + 1;
% end
% C = (169).*ones(size(node,1),3);
% %C(:,1) = 0.671;
% %C(:,2) = 0.329;
% C(pial_conflict,1) = 255;
% C(pial_conflict,2) = 0;
% C(pial_conflict,3) = 0;
% writeOFF('lh_pial_noCC_color.off', node,face,[],C,[]);
% img=figure(2);
%     clf 
%     set(img, 'Position', [300 100 600 600]); 
%     hold on
%     %axis equal
%     set(gca, 'Layer', 'top')
%     pp= patch('Faces',face,'Vertices',node);
%     view(180, -90);
%     %camroll(90)
%     ax = gca;
%     %ax.CameraPosition = [0 10 0];
%     caxis([0,400]);
%     axis off;
% %     %colormap jet
% %     %colorbar
%     hold on
%     sphereSize = 0.008; 
%     c = [1 0 0];
%     drawSpheres(node(Lmks,:),sphereSize,c);
    
save(fullfile(base_path,dire(i).name,'fan_landmark1000.txt'),'Lmks','-ascii');
end
%end
%end
