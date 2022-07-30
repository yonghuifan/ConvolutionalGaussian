function drawSpheres(X,r,c)

if size(c,1)==1
    c = repmat(c,[size(X,1),1]);
end
N = 30;
[x,y,z] = sphere(N);
% S = repmat(100,[size(X,1),1]);
% C = repmat(3,[size(X,1),1]);
% s = S(:);
% c = C(:);
% r = r + 0.1;
for ii = 1:size(X,1)
    surf(r*x+X(ii,1), r*y+X(ii,2), r*z+X(ii,3),'edgecolor','none','facecolor',c(ii,:));
    %scatter3(X(ii,1), X(ii,2), X(ii,3),180,'LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','r');
    %light('Position',[0, 0, X(ii,3)],'Style','local');
    %patch(r*x+X(ii,1), r*y+X(ii,2), r*z+X(ii,3),'FaceLighting','phong','edgecolor','none','facecolor',c(ii,:));
    %camlight
    %light('Position',[X(ii,1), X(ii,2), X(ii,3)]);
    %lighting phong
    %camlight
    %hold on
end
end