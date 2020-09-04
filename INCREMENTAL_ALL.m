function [S, time] = INCREMENTAL_ALL(A, node_label1, edge_label1, Q, node_label2, edge_label2, ...
    Qr, node_label3, edge_label3, H, l, alpha, l_node, l_edge, r, t)
%    FIRST-E Implementation:
%   If only edge attribute is changed between query revision, l is not
%   empty and serves as the index of revised edges, otherwise l is an empty
%   vector; The attributed network: A, node_label1, edge_label1; The
%   initial query network: Q, node_label2, edge_label2; The revised query
%   network: Qr, node_label3, edge_label3;
%   l_node: number of node attribute;
%   l_edge: number of edge attribute;
%   r: top r eigen decomposition of element-wise product of 'EA' and A;
%   t: top t eigen decomposition of element-wise product of 'EQ' and Q;

%%  Precomputing Stage:

m = size(A,1);
n = size(Q,1);
H=H./sum(sum(H)); h = H(:); LAM =[];
UA = cell(1, l_edge); 
for i = 1:l_edge
    UA{i} = zeros(n,r);
end
LamA = cell(1, l_edge); 
for i = 1:l_edge
    LamA{i} = zeros(r,r);
end
Uq = cell(1, l_edge); 
for i = 1:l_edge
    Uq{i} = zeros(m,t);
end
Lamq = cell(1, l_edge); 
for i = 1:l_edge
    Lamq{i}=zeros(t,t);
end
UQ = cell(1, l_edge); 
for i = 1:l_edge
    UQ{i} = zeros(m,t);
end
LamQ = cell(1, l_edge); 
for i = 1:l_edge
    LamQ{i} = zeros(t,t);
end
deg = []; E1 = cell(l_edge, 1); E2 = cell(l_edge, 1); N = [];
U = zeros(m*n, l_edge*t*r);
for i = 1 : l_edge
    [re1, ce1, ~] = find(edge_label1 == i);
    E1{i} = sparse(re1, ce1, 1, m, m);
    [UA{i}, LamA{i}] = eigs(E1{i}.*A, r);
end
for i = 1:l_edge
    [re1, ce1, ~] = find(edge_label2 == i);
    E2{i} = sparse(re1, ce1, 1, n, n);
    [Uq{i}, Lamq{i}] = eigs(E2{i}.*Q, t);
end
E2 = cell(l_edge, 1);
%%  Interactive Stage:
tic;
if isempty(l)
% Construct N and D;    
    for i = 1 : l_node
        [rn1, cn1, ~] = find(node_label1 == i);
        [rn2, cn2, ~] = find(node_label3 == i);
        N1 = sparse(rn1, cn1, 1, m, 1);
        N2 = sparse(rn2, cn2, 1, n, 1);
    if isempty(N), N = kron(N2, N1);
    else N = N + kron(N2, N1); end
    end
    
    for i = 1:l_edge
        [re2, ce2, ~] = find(edge_label3 == i);
        E2{i} = sparse(re2, ce2, 1, n, n);
        for j = 1:l_node
            [rn1, cn1, ~] = find(node_label1 == j);
            [rn2, cn2, ~] = find(node_label3 == j);
            N1 = sparse(rn1, cn1, 1, m, 1);
            N2 = sparse(rn2, cn2, 1, n, 1);
            if isempty(deg), deg = kron((Qr.*E2{i})*N2,(A.*E1{i})*N1);
            else deg = deg + kron((Qr.*E2{i})*N2,(A.*E1{i})*N1); end
        end
    end
% approximate element-wise product of EQ and A;
    for j = 1:l_edge
        [UQ{j}, LamQ{j}] = eigs(E2{j}.*Qr, t);
    end
    for i = 1:l_edge
        LAM = blkdiag(LAM, kron(LamQ{i},LamA{i}));
        U(1:m*n, (t*r*(i-1)+1):(t*r*i)) = kron(UQ{i}, UA{i});
    end
    clear UA;
    clear UQ;
    deg2 = deg.^(0.5);
    deg2(deg2==Inf) = 0;
    deg2 = full(deg2);
    N = full(N);
    S = bsxsolver_mex(deg2, N, U, LAM, alpha, h, m, n);
%      L = bsxfun(@times, deg2, (bsxfun(@times, N, U)));
%      R = bsxfun(@times, bsxfun(@times, U', N'), deg2');
%     clear U;
%     iLAM = pinv(LAM);
%     temp = pinv(iLAM-alpha*R*L);
%     s = (1-alpha)*(h+alpha*(L*(temp*(R*h))));
%     S = reshape(s,m,n);
else
    % Construct N and D;    
    for i = 1 : l_node
        [rn1, cn1, ~] = find(node_label1 == i);
        [rn2, cn2, ~] = find(node_label2 == i);
        N1 = sparse(rn1, cn1, 1, m, 1);
        N2 = sparse(rn2, cn2, 1, n, 1);
    if isempty(N), N = kron(N2, N1);
    else N = N + kron(N2, N1); end
    end
    
    for i = 1:l_edge
        [re2, ce2, ~] = find(edge_label3 == i);
        E2{i} = sparse(re2, ce2, 1, n, n);
        for j = 1:l_node
            [rn1, cn1, ~] = find(node_label1 == j);
            [rn2, cn2, ~] = find(node_label2 == j);
            N1 = sparse(rn1, cn1, 1, m, 1);
            N2 = sparse(rn2, cn2, 1, n, 1);
            if isempty(deg), deg = kron((Q.*E2{i})*N2,(A.*E1{i})*N1);
            else deg = deg + kron((Q.*E2{i})*N2,(A.*E1{i})*N1); end
        end
    end
% only need to approximate element-wise product of EQ and A in the index l;
    In = [];
    for j = size(l,1)
        In = [In, edge_label2(l(j,1),l(j,2)), l(j, 3)];
    end
    In = unique(In);
    for j = 1:size(In,1)
        [Uq{In(j)}, Lamq{In(j)}] = eigs(E2{In(j)}.*Qr, t);
    end
    for i = 1:l_edge
        LAM = blkdiag(LAM, kron(Lamq{i},LamA{i}));
        U(1:m*n, (t*r*(i-1)+1):(t*r*i)) = kron(Uq{i}, UA{i});
    end
    clear UA;
    clear Uq;
    deg2 = deg.^(0.5);
    deg2(deg2==Inf) = 0;
    deg2 = full(deg2);
    N = full(N);
    S = bsxsolver_mex(deg2, N, U, LAM, alpha, h, m, n);
%      L = bsxfun(@times, deg2, (bsxfun(@times, N, U)));
%      R = bsxfun(@times, bsxfun(@times, U', N'), deg2');
%     clear U;
%     iLAM = pinv(LAM);
%     temp = pinv(iLAM-alpha*R*L);
%     s = (1-alpha)*(h+alpha*(L*(temp*(R*h))));
%     S = reshape(s,m,n);
end
time = toc;
fprintf('running time = %f\n', time);

end