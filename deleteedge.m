function G = deleteedge(G, paths)
l = size(paths,2);
for i = 1:l-1
    edge = [];
%    edge = cell2mat(paths);
    G(paths(i), paths(i+1)) = 0;
    G(paths(i+1), paths(i)) = 0;
end
end