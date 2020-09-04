function [SS,time] = INCREMENTAL_E(A, Q, H, node_label1, node_label2, l_node, edge_label1, ...
edge_label2, l_edge, alpha, r, relax, l, a, b, firstrun)
%% Consider edge attribute but only edge attribute changes in the incremental process
%  l, (j,k) are index of changed edge attribute and index of the edge
%  with changed attribute respectively
%  If m, n are both huge, Bin may require huge memory spaces that exceed
%  maximum!
m = size(A,1);
n = size(Q,1);
% [UA,LambdaA]=eigs(A,r);
% [UQ,LambdaQ]=eigs(Q);
H=H./sum(sum(H));y = H(:);
% compute D
deg = []; E1 = cell(l_edge, 1); E2 = cell(l_edge, 1); E=[]; N = [];
for i = 1:l_node
    [rn1, cn1, ~] = find(node_label1 == i);
    [rn2, cn2, ~] = find(node_label2 == i);
    N1 = sparse(rn1, cn1, 1, m, 1);
    N2 = sparse(rn2, cn2, 1, n, 1);
    if isempty(N), N = kron(N2, N1);
    else N = N + kron(N2, N1); end
end

for i=1:l_edge
    [re1, ce1, ~] = find(edge_label1 == i);
    [re2, ce2, ~] = find(edge_label2 == i);
    E1{i} = sparse(re1, ce1, 1, m, m);
    E2{i} = sparse(re2, ce2, 1, n, n);
    if relax == 0
        for j = 1:l_node
            [rn1, cn1, ~] = find(node_label1 == j);
            [rn2, cn2, ~] = find(node_label2 == j);
            N1 = sparse(rn1, cn1, 1, m, 1);
            N2 = sparse(rn2, cn2, 1, n, 1);
            if isempty(deg), deg = kron((Q.*E2{i})*N2,(A.*E1{i})*N1);
            else deg = deg + kron((Q.*E2{i})*N2,(A.*E1{i})*N1); end
        end
    end
    if isempty(E)
        E = kron(E2{i},E1{i});
    else
        E = E + kron(E2{i},E1{i});
    end
end
if relax == 1, deg = kron(sum(Q,2),sum(A,2)); end
    
D = 1./sqrt(deg);
D(D == Inf) = 0;

% computing U, V, Lambda, B and all unchanged variables
[U1,Lambda1]=eigs(E1{l},r);
X = bsxfun(@times, D, N);
X = N;
if firstrun == 1    
   I = sparse(1:m*n,1:m*n,1,m*n,m*n);
   B = I-alpha.*bsxfun(@times, X, bsxfun(@times, (E.*(kron(Q,A))),X')); 
   Bin = inv(B);% takes lots of time, but can be pre-computed (firstrun=0)
   Bin(Bin==Inf)=0;
%     EAQ = bsxfun(@times, X', bsxfun(@times, (E.*(kron(A,Q))),X));
%     
%     [UB, LambdaB] = eigs(EAQ,r);
%     
%     In = inv(inv(LambdaB) - alpha.*UB'*UB);
%     Bin = I + alpha.*UB*In*UB';
    save('Bin(DBLP).mat', 'Bin');
else
    load('Bin(DBLP).mat');
end
tic;
if firstrun == 1
    s = (1-alpha).*(Bin*y);
    SS = reshape(s,m,n);
    fprintf('running time = %f\n',toc);
    time = toc;
else
    Eq1 = sparse(a,b,1,n,n);
    Eq1 = Eq1 + Eq1';
    [U2,Lambda2]=eigs(Eq1);
    U = bsxfun(@times, X, kron(U2,U1));
    V = bsxfun(@times, kron(U2',U1'), X');
    Lambda = kron(Lambda2,Lambda1);

    % computing new similarity matrix
    eta = pinv(Lambda)-alpha.*V*Bin*U;
    eta = pinv(eta);
    s = (1-alpha).*(Bin*y+alpha.*Bin*(U*(eta*(V*(Bin*y)))));
    SS = reshape(s,m,n);
    fprintf('running time = %f\n',toc);
    time = toc;
end
end