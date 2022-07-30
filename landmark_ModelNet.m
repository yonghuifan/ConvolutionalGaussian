%%% init
clear
close all;
initialize;
yalmip('clear');
rng(1);

fid = fopen('G:\modelnet40_normal_resampled\piano\piano_0005.txt');
% fid = fopen('G:\modelnet40_normal_resampled\airplane\airplane_0008.txt');
% fid = fopen('G:\SHREC14\Real\Data\0.off');
[A,cnt] = fscanf(fid,'%f,%f,%f,%f,%f,%f\n');
A = reshape(A, 6,[]);
node = A(1:3,:);
normal = A(4:6,:);
fclose(fid);
pcshow(node');
G.V = node;
G.nV = 10000;
params.preprocess.dataFileType = 'off';
numLmk = 500;
BNN = 100;

for i = 1:15

    %dist = load(fullfile(base_path,[dire(i).name,'_geodesicDist.txt']));
    %idx = load(fullfile(base_path,[dire(i).name,'_geodesicIdx.txt']));
    [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
    color = zeros(10000,3);
    % color(Lmks,1) = 100;
    shape = ShapeStruct(G.V', []);
    %[ fig, S ] = RenderPointCloudFunction(shape, color);    
    Lmks = reshape(Lmks,1,[]);
     
    %save([dataFolder,dataFileNames{j},'_landmark.txt'],'tmp','-ascii');
    
    % [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
    node = G.V';
    face = G.F';
    %[BV,BE] = FindBoundaries(G);
    %tmp = setdiff(sort(Lmks), sort(BV));
    tmp = Lmks;
    
    CMap=colormap(parula(10));
    [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
    %[GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
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
    pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
    %pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
    %view(90,180);
    ax = gca;
    %ax.CameraPosition = [0 10 0];
    caxis([0,400]);
    axis off;
    %colormap jet
    %colorbar
    hold on
    sphereSize = 0.004; 
    c = [1 0 0];
    drawSpheres(node(tmp,:),sphereSize,c);
    % drawSpheres(node(Lmks{j},:),sphereSize,[1 0 1]);
end
