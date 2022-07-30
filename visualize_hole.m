%%% init
clear
close all;
initialize;
yalmip('clear');
rng(1);

base_path = 'G:\CPsurfcomp\DATA\teeth\meshes\w09_sas.off';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\mt1\meshes';
%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\radius\meshes';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC14\Real\Data';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\GPLmkBDMatch-master\sample_data';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

params.preprocess.dataFileType = 'off';

G = Mesh(params.preprocess.dataFileType,fullfile(base_path));
node = G.V;
node = PerformMeshSmoothing(G,node);
G.V = node';
G.Normalize();
[~,~,flip] = G.ComputeNormal();
if flip
    G.F = G.F([1 3 2],:);
end
% [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
node = G.V';
face = G.F';
[BV,BE] = FindBoundaries(G);
%[C, ia, ib] = intersect()
%tmp = setdiff(Lmks, sort(BV),'stable');
%tmp(Lmk_needed+1:end) = [];
% tmp = Lmks;
%save([fullfile(base_path,[dire(i).name]),'_SM_landmark1000.txt'],'tmp','-ascii');
%CMap=colormap(parula(10));
% [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
%[GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
%MC = MC;

% MC_min = min(MC);
% MC_max = max(MC);
% GC_min = min(GC);
% GC_max = max(GC);
node(:,3) = -node(:,3);
%img=figure(1);
%clf
%set(img, 'Position', [300 100 600 600]);
%hold on
%axis equal
%set(gca, 'Layer', 'top')
%pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(GC*2),'FaceColor','interp','EdgeColor','none');
%pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
%view(90,180);
%ax = gca;
%ax.CameraPosition = [0 10 0];
%caxis([0,400]);
%axis off;
%colormap jet
%colorbar
%hold on
%sphereSize = 0.01;
c = [1 0 0];
%drawSpheres(node(BV,:),sphereSize,c);
color = ones(size(node,1),3);
color(:,1) = 227;
color(:,2) = 218;
color(:,3) = 201;
color(BV,1) = 255;
color(BV,2) = 0;
color(BV,3) = 0;
writeOFF('bone_vis.off', node, face, [], color);


