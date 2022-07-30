function [node,face] = read_oasis(filename)

fid = fopen(filename);

A = fscanf(fid, 'Vertex %d %f %f %f {wid=%d normal=(%f %f %f)}\n');
A = reshape(A, [8, length(A)/8]);
vertex = A(2:4,:);
normal = A(6:8,:);
index = A(1,:);
% read faces
A = fscanf(fid, 'Face %d %d %d %d {matid=%d}\n');
A = reshape(A, [4, length(A)/4]);
face = A(2:4,:);

for i = 1:size(index,1)
    [row,col] = find(face == index(i));
       for m = 1:size(row,1)
            face(row(m),col(m)) = i;
       end
end
fclose(fid);

end

