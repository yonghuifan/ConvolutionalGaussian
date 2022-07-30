clear;close all;

start_path = 'G:\OASIS\pos';

for n = 1:52
    tic
    left_name = fullfill(start_path, 'l', [int2str(n),'.txt']);
    right_name = fullfill(start_path, 'r', [int2str(n),'.txt']);
    [lnode,lface] = read_oasis(left_name);
    [rnode,relem]=read_oasis(right_name);
    
    output_name = fullfill(start_path, [int2str(n),'.txt']);
    
    toc
end