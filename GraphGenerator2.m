function [EDGES, NODES, MatchNode] = GraphGenerator2(S,A,Q,l)
%% Generate Subgraph directly from S using top p exhaustive method 
% -- INPUTS: 
%         S: The cross-network similarity between query and data graph
%         A: The adjacency matrix of data network
%         Q: The adjacency matrix of query
%         l: The size of searching space
% 
% -- OUTPUTS:
%         EDGES: edge list of the matching subgraph
%         NODES: node list of the matching subgraph
%         MatchNode: matching nodes in the data network
% calculate X (k*m 0,1 matrix)

if size(S,1)>size(S,2)
    S=S';
end
m = size(S,1);
n = size(S,2);
a = 0.5;
b = 0.5;
I = diag(ones(1, m));
[~, index] = sort(S,2,'descend');
index_select = index(:,1:l);
perm_result = [];
total = l.^m;
for i = 0:total - 1
    val = i;
    info = zeros(1,m);
    for j = m:-1:1
        info(j) = mod(val,l);
        val = fix(val./l);
    end
    tem = [];
    for k = 0:m-1
        tem = [tem, index_select(k+1, info(k+1)+1)];
    end
    perm_result = [perm_result; tem];
end

FX = [];
for i = 1:size(perm_result,1)
    test_x = zeros(m,n);
    for j = 1:m
        test_x(j,perm_result(i,j)) = 1;
    end
    test_x = sparse(test_x);
    [fx, ~]=GoodnessF(test_x,A,Q,I,S,a,b);
    FX = [FX, fx];
end
[~,i] = min(FX);
final_x_index = perm_result(i,:);
X = sparse(1:m, final_x_index, ones(1,m), m, n);

% Generate nodes and edges
Q=sparse(Q);
T = X'*Q*X - and(A, X'*Q*X);
[row,col] = find(tril(T) == 1);
paths = cell(size(row,1));
[row3, col3] = find(tril(and(A, X'*Q*X)));
% nodes = [(1:size(A,1)); 100*rand(2, size(A,1))]';
% num_edge = size(find(A == 1))/2;
% [row2,col2] = find(tril(A) == 1);
% segments = [(1:num_edge); row2'; col2']';
% for i = 1:size(row)/2
%     [~, path] = dijkstra(nodes, segments, row(i,1), col(i,1));
%     if isempty(shortest_p)
%         shortest_p = path;
%     else
%         shortest_p = {shortest_p path};
%         
%     end
% end
Z = A;
Z(row3, col3) = 0;
Z(col3, row3) = 0;
for i = 1:size(row,1)
        [~,path] = Dijk4(Z, Z, row(i,1), col(i,1));
        paths(i,i) = {path};
        Z=deleteedge(Z,path);
end
clear Z;
    
%[row3, col3] = find(tril(and(A, X'*Q*X)));
select_p = [];
if size(paths,1)==1
    select_p = cell2mat(paths);
else
    for i = 1:size(paths,1)
        select_p = cat(2, select_p, paths(i,i));
    end
    select_p = cell2mat(select_p);
end

nodes = cat(2, select_p, row3', col3');
% nodes = cell2mat(nodes);
NODES = unique(nodes);
e = [];

for i = 1:size(paths,1)
    if size(paths,1)==1
        p = cell2mat(paths);
    else
        p = cell2mat(paths(i,i));
    end
    for j = 1:size(p,2)-1
        if isempty(e)
            e = [p(1,j),p(1,j+1)];
        else
            e = [e; [p(1,j), p(1,j+1)]];
        end
    end
end
new_edge = [[row3, col3]; e];
EDGES = new_edge;
MatchNode = final_x_index;
end