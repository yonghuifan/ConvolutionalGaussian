%%% init
clear
close all;
%initialize;
% yalmip('clear');
rng(1);

%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\teeth\meshes';
base_path = 'G:\ADNI2_processed\AD\002_S_5018\';
[pv,pf]=readoff([base_path,'\tpial_fixed.off']);
[wv,wf]=readoff([base_path,'\twhite_fixed.off']);
% intersection detection and mark
[P, w_index, p_index] = intersect(wv,pv,'rows');    
[pn,connnum,count]=meshconn(pf,size(pv,1));
[wn,connnum,count]=meshconn(wf,size(wv,1));
    
        count = 1;
        remove_pv = p_index';
        remove_wv = w_index';
        while count < 3
            for i = 1:size(p_index,1)
                remove_pv = [remove_pv, unique(pn{remove_pv(1,i)})];
            end
            remove_pv = unique(remove_pv);
            count = count + 1;
        end
    
        while count < 6
            for i = 1:size(w_index,1)
                remove_wv = [remove_wv, unique(wn{remove_wv(1,i)})];
            end
            remove_wv = unique(remove_wv);
            count = count + 1;
        end

dire = dir(fullfile(base_path, '*.off'));

params.preprocess.dataFileType = 'off';
numLmk = 20;

BNN = 1800;

% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

for i = 7:length(dire)% 7 28 26
    
    %if ~exist(outputfilename,'file')
    display(dire(i).name);
    G = Mesh(params.preprocess.dataFileType,fullfile(base_path,[dire(i).name]));
    node = G.V;
    %node = PerformMeshSmoothing(G,node);
    %G.V = node';
    G.Normalize();
    [~,~,flip] = G.ComputeNormal();
    if flip
        G.F = G.F([1 3 2],:);
    end
    
    %dist = load(fullfile(base_path,[dire(i).name,'_geodesicDist.txt']));
    %idx = load(fullfile(base_path,[dire(i).name,'_geodesicIdx.txt']));
    %dist = dist(:,1:BNN+1);
    %idx = idx(:, 1:BNN+1);
    [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
    %Lmks = GetGPLmk(G, numLmk);
    % Lmks = GetGPLmk1(G, numLmk, BNN, idx, dist);
    % Lmks = GetGPLmk_Euclidean(G, numLmk);
    %Lmks = GetGPLmk_NoCurvature(G, numLmk);
    % Lmks = GetGPLmk_periodic(G, numLmk, BNN, idx, dist);
    Lmks = reshape(Lmks,1,[]);
     
    % [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
    node = G.V';
    face = G.F';
    %[BV,BE] = FindBoundaries(G);
    %[C, ia, ib] = intersect()
    %tmp = setdiff(Lmks, sort(BV),'stable');
    %tmp(Lmk_needed+1:end) = [];
    tmp = Lmks;
    %parsave([fullfile(base_path,[dire(i).name]),'_periodic_landmark1000.txt'],tmp);
    CMap=colormap(parula(10));
    % [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
    [GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
    %MC = MC;
    GC = normalize(GC);
    % MC_min = min(MC);
    % MC_max = max(MC);
    % GC_min = min(GC);
    % GC_max = max(GC);
    node(:,3) = -node(:,3);
    img=figure(1);
    clf 
    set(img, 'Position', [300 400 400 400]); 
    hold on
    %axis equal
    set(gca, 'Layer', 'top')
    pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(GC*400),'FaceColor','interp','EdgeColor','none','BackFaceLighting','lit');
    %pp= patch('Faces',face,'Vertices',node,'FaceColor','interp','EdgeColor','none');
    view(-132, 44);
    ax = gca;
    %ax.CameraPosition = [0 10 0];
    caxis([0,400]);
    axis off;
    %colormap jet
    %colorbar
    hold on
    sphereSize = 0.018; 
    c = [1 0 0];
    
    drawSpheres(node(tmp,:),sphereSize,c);
%     color = ones(size(node,1),3);
%     color(:,1) = 227;
%     color(:,2) = 218;
%     color(:,3) = 201;
%     color(BV,1) = 255;
%     color(BV,2) = 0;
%     color(BV,3) = 0;
    %writeOFF('bone_vis.off', node, face, [], color);
    % end
end
