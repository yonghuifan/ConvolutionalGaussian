clear
close all;

base_path = 'G:\OASIS\neg\l';

dire = dir(base_path);
dire(1:2) = [];

for i = 1:length(dire)
    if strfind(dire(i).name,'.off') 
        filename = fullfile(dire(i).folder, dire(i).name);
        [node,face] = readoff(filename);
        DT = triangulation(face, node);
        normal = vertexNormal(DT);
        node = [node, normal];
        outfilename = fullfile(dire(i).folder, [int2str(i),'.txt']);
        save(outfilename, 'node', '-ascii', '-double', '-tabs');
    end
%     fileID = fopen(outfilename,'w');
%     fprintf(fileID,'%f %f %f %f %f %f\n', node);
%     fclose(fileID);
end