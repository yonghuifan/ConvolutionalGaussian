function oasis_m_to_off(filename, output_name)
fid = fopen(filename);

A = fscanf(fid, 'Vertex %d %f %f %f {Jfeature=(%f %f %f %f %f %f %f)}\n');
A = reshape(A, [11, length(A)/11]);
vertex = A(2:4,:);
%uv = A(5:6,:);
% read faces
A = fscanf(fid, 'Face %d %d %d %d\n');
A = reshape(A, [4, length(A)/4]);
face = A(2:4,:);

fclose(fid);
writeOFF(output_name,vertex',face');
end