function [distBD, Lmks, putativeInds, matchIndsBD, regParamBD, isoParamBD] = wrapperLmkMatch(Gs,params,dataFolder,dataFileNames)

% compute GP landmarks
for j = 1:2
    switch params.typeLmk
        case 'GP'
            BNN = 100;
            %                         if j == 1
            %                         dist = load('C:\Users\yfan61\Documents\MATLAB\GPLmk1\sample_data\a16_sas_aligned_geodesicDist.txt');
            %                         idx = load('C:\Users\yfan61\Documents\MATLAB\GPLmk1\sample_data\a16_sas_aligned_geodesicIdx.txt');
            %                         dist = dist(:,BNN+1);
            %                             idx = idx(:,BNN+1);
            %                         else
            %                             dist = load('C:\Users\yfan61\Documents\MATLAB\GPLmk1\sample_data\b02_sas_aligned_geodesicDist.txt');
            %                             idx = load('C:\Users\yfan61\Documents\MATLAB\GPLmk1\sample_data\b02_sas_aligned_geodesicIdx.txt');
            %                             dist = dist(:,BNN+1);
            %                             idx = idx(:,BNN+1);
            %                         end
            %             display('Start calculating geodesic distance...');
            %             for i = 1:Gs{j}.nV
            %                 [D,S,Q] = PerformFastMarching(Gs{j}, i);
            %                 [dist(i,:), idx(i,:)] = mink(D, BNN+1);
            %             end
            [idx, dist] = knnsearch(Gs{j}.V',Gs{j}.V','K',BNN+1);
            Lmks{j} = Gs{j}.GetGPLmk1(params.GPLmk.numLmk, BNN, idx, dist);
            %Lmks{j} = Gs{j}.GetGPLmk_Fan(params.GPLmk.numLmk, BNN, idx, dist);
        case 'GP_nW'
            Lmks{j} = Gs{j}.GetGPLmk_NoCurvature(params.GPLmk.numLmk);
        case 'GP_Euc'
            Lmks{j} = Gs{j}.GetGPLmk_Euclidean(params.GPLmk.numLmk);
        case 'GT'
            Lmks{j} = Gs{j}.Aux.GTLmks(1:params.GTLmk.numLmk);
        otherwise
            error('invalid typeLmk')
    end
    Lmks{j} = reshape(Lmks{j},1,[]);
    Lmks{j} = Lmks{j}(160:end);
    %tmp = Lmks{j};
    %     o = full(outline_loop((Gs{j}.F)'));
    %     tmp = setdiff(Lmks{j},o);
    %save([dataFolder,'\',dataFileNames{j},'_LM2.txt'],'tmp','-ascii');
end
%for j = 1:2
%    outputfile =[dataFolder,'\',dataFileNames{j},'_V2.off'];
%    virsualize_surface(Lmks{j},Gs{j}.V',Gs{j}.F',outputfile);
%end
% compute putative correspondences
[putativeMatches, putativeInds] = computePutativeMatches(Gs,Lmks,params.computePutativeMatches);
% compute BD filtered landmark matches
[matchIndsBD,regParamBD,isoParamBD] = matchSurfaceLmksBD(Gs,Lmks,putativeMatches,params.matchSurfaceLmksBD);
% compute the Induced CP Distance
distBD = MapToDist(Gs{1}.V,Gs{2}.V,...
    knnsearch(regParamBD{2}.V(1:2,:)',regParamBD{1}.V(1:2,:)'),...
    Gs{1}.Aux.VertArea);
fprintf('Lmk-BD: Induced CP distance %f\n', distBD);