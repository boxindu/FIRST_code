function [S, time] = INCREMENTAL_N(A, Q, H, node_label1, node_label2, node_labelq, l_node, alpha,r, relax, issecondrun)
%% Consider node attribute but only N changes in the incremental process
%  node_labelq is the node attribute vector after node attribute change;
%  issecondrun takes 0 or 1, indicating wether the current run is before or
%  after node attribute change.
m = size(A,1);
n = size(Q,1);
% Low rank approximation on two input adjacency matrix A, Q
[UA,LambdaA]=eigs(A,r);
[UQ,LambdaQ]=eigs(Q,5);
% compute D 
H=H./sum(sum(H));y = H(:);
N = []; deg = [];
if issecondrun == 0
for i = 1:l_node
    [rn1,cn1,~] = find(node_label1 == i);
    [rn2,cn2,~] = find(node_label2 == i);
    N1 = sparse(rn1, cn1, 1, m, 1);
    N2 = sparse(rn2, cn2, 1, n, 1);
    if isempty(N), N = kron(N2, N1);
    else N = N + kron(N2, N1); end 
    if relax == 0
        if isempty(deg), deg = kron(Q*N2, A*N1);
        else deg = deg + kron(Q*N2, A*N1); end
    end
end
if relax == 1, deg = kron(sum(Q,2),sum(A,2)); end

else
    tic;
    for i = 1:l_node
    [rn1,cn1,~] = find(node_label1 == i);
    [rn2,cn2,~] = find(node_labelq == i);
    N1 = sparse(rn1, cn1, 1, m, 1);
    N2 = sparse(rn2, cn2, 1, n, 1);
    if isempty(N), N = kron(N2, N1);
    else N = N + kron(N2, N1); end 
    if relax == 0
        if isempty(deg), deg = kron(Q*N2, A*N1);
        else deg = deg + kron(Q*N2, A*N1); end
    end
    end
if relax == 1, deg = kron(sum(Q,2),sum(A,2)); end
    t1 = toc;
end

tic;
D=1./sqrt(deg);
D(D==Inf)=0;
X = D.*N;
iD2 = X.^2;
iX = 1./X;
iX(iX==Inf)=0;

U2 = kron(UQ, UA);% Large memory space needed, sparse doesn't work(too slow)
clear UA;
C2 = kron(LambdaQ, LambdaA);
clear N;
clear N1;
clear rn1;
clear cn1;
clear X;
clear deg;
clear D;
clear H;
eta = U2'*bsxfun(@times,iD2,U2);
temp = sparse(pinv(pinv(C2)-alpha.*eta));% 
clear eta;
clear C2;
iXD2 = iX.*iD2;
s = (1-alpha).*(y+alpha.*(iXD2.*(U2*(temp*(U2'*(iXD2.*y))))));
clear temp;
clear iD2;
clear U2;
t2 = toc;
fprintf('running time = %f\n', t1+t2);
time = t1+t2;
S=reshape(s,m,n);
end