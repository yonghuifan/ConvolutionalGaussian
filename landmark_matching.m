clear;close all;
dataFolder = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\mt1\meshes\';
dire = dir(fullfile(dataFolder''));
params.GTLmk.numLmk = 16;
dataFileNames = {'Aa211457_sas', 'Aa215050_sas'};
params.preprocess.force = false;
params.preprocess.suffix = '_PRE';
params.preprocess.dataFileType = 'off';
params.preprocess.CPM.NumDensityPnts = 100;

params.CPM.FeatureType = 'ConfMax';
params.CPM.AngleIncrement = 0.05;
params.CPM.NumFeatureMatch = 4;
params.CPM.GaussMinMatch = 'off';

params.GPLmk.numLmk = 40;
params.computePutativeMatches.forceIdentity = false; % enabling makes sense only for GT landmarks
params.computePutativeMatches.WKSNN = 2;
%params.matchSurfaceLmksBD.visualize = visualize;
params.matchSurfaceLmksBD.forceIdentity = false; % enabling makes sense only for GT landmarks
params.matchSurfaceLmksBD.paramBDfilt.K = 1.5;
params.matchSurfaceLmksBD.paramBDfilt.pnorm = 0.1;
params.matchSurfaceLmksBD.paramBDfilt.initdelta = 1;
params.matchSurfaceLmksBD.paramBDfilt.boundaryConstraintType = 'similarity+translation'; %{'none','fixed','linear','similarity','affine','similarity+translation'}

params.visualizeCorrepsondences2D.offsetFactors = [-1.1 0];
params.visualizeCorrepsondences2D.rotationAngle = 0*(pi/180);

% update parameters according to selected method
for i = 1:2
    Lmks{i} = load([dataFolder, dataFileNames{i},'.off_Fan2_landmark1000.txt']);
    Lmks{i} = Lmks{i}(1:100);
    CMap=colormap(parula(10));
    G{i} = Mesh(params.preprocess.dataFileType,fullfile(dataFolder,[dire(i).name]));
    G{i}.Normalize();
    % reorient (outward facing normals)
    [~,~,flip] = G{i}.ComputeNormal();
    if flip
        G{i}.F = G{i}.F([1 3 2],:);
    end
    % compute WKS
    G{i}.Aux.WKS = G{i}.ComputeWKS([]);
    [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G{i});
    node = G{i}.V';
    face = G{i}.F';
    [GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
    %MC = MC;
    % MC_min = min(MC);
    % MC_max = max(MC);
    % GC_min = min(GC);
    % GC_max = max(GC);
    node(:,3) = -node(:,3);
    img=figure(i);
    clf 
    set(img, 'Position', [300 100 600 600]); 
    hold on
    %axis equal
    set(gca, 'Layer', 'top')
    pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(GC),'FaceColor','interp','EdgeColor','none');
    %pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
    %view(90,180);
    ax = gca;
    %ax.CameraPosition = [0 10 0];
    caxis([0,400]);
    axis off;
    %colormap jet
    %colorbar
    hold on
    sphereSize = 0.007; 
    c = [1 0 0];
    drawSpheres(node(Lmks{i},:),sphereSize,c);
end

[putativeMatches, putativeInds] = computePutativeMatches(G,Lmks,params.computePutativeMatches);
[matchIndsBD,regParamBD,isoParamBD] = matchSurfaceLmksBD(G,Lmks,putativeMatches,params.matchSurfaceLmksBD);
distBD = MapToDist(G{1}.V,G{2}.V, knnsearch(regParamBD{2}.V(1:2,:)',regParamBD{1}.V(1:2,:)'),G{1}.Aux.VertArea);
fprintf('Lmk-BD: Induced CP distance %f\n', distBD);
[distCPM, regParamCP, CPM12] = wrapperCPM(G,params);






