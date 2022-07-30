%%% init
clear
close all;
initialize;
yalmip('clear');
rng(1);

%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\teeth\meshes';
base_path = 'G:\Barn.ply';
%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\radius\meshes';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC14\Real\Data';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\GPLmkBDMatch-master\sample_data';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

params.preprocess.dataFileType = 'ply';
numLmk = 50;
Lmk_needed = 50;
BNN = 1000;
%outputfilename = [fullfile(base_path,[dire(i).name]),'_Fan_landmark1000.txt'];
%if ~exist(outputfilename,'file')
%display(dire(i).name);
%G = Mesh(params.preprocess.dataFileType,base_path);
%node = G.V;
%node = PerformMeshSmoothing(G,node);
ptCloud = pcread(base_path);
G.V = ptCloud.Color;
% G.Normalize();
% [~,~,flip] = G.ComputeNormal();
% if flip
%     G.F = G.F([1 3 2],:);
% end

%     dist = load(fullfile(base_path,[dire(i).name,'_geodesicDist.txt']));
%     idx = load(fullfile(base_path,[dire(i).name,'_geodesicIdx.txt']));
%     dist = dist(:,1:BNN+1);
%     idx = idx(:, 1:BNN+1);
[idx, dist] = knnsearch(G.V,G.V,'K',BNN+1);
Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
%     % Lmks = GetGPLmk(G, numLmk);
%     % Lmks = GetGPLmk_Euclidean(G, numLmk);
%     % Lmks = GetGPLmk_NoCurvature(G, numLmk);
%
%     Lmks = reshape(Lmks,1,[]);
%
%     % [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
node = G.V';
face = G.F';
%     [BV,BE] = FindBoundaries(G);
%     %[C, ia, ib] = intersect()
%     tmp = setdiff(Lmks, sort(BV),'stable');
%     tmp(Lmk_needed+1:end) = [];
tmp = Lmks;
%save([fullfile(base_path,[dire(i).name]),'_Fan_landmark1000.txt'],'tmp','-ascii');
CMap=colormap(parula(10));
[Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
% [GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
%MC = MC;

% MC_min = min(MC);
% MC_max = max(MC);
% GC_min = min(GC);
% GC_max = max(GC);
node(:,3) = -node(:,3);
img=figure(1);
clf
set(img, 'Position', [300 100 600 600]);
hold on
%axis equal
set(gca, 'Layer', 'top')
pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(Cgauss*2),'FaceColor','interp','EdgeColor','none');
%pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
%view(90,180);
ax = gca;
%ax.CameraPosition = [0 10 0];
caxis([0,400]);
axis off;
%colormap jet
%colorbar
hold on
sphereSize = 0.006;
c = [1 0 0];
drawSpheres(node(tmp,:),sphereSize,c);
% end
