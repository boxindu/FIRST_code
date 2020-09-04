function [SS, time] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,r, relax, issecondrun)
%% Consider node attribute and do not consider edge attribute, but only Q changes in the incremental process
m = size(A,1);
n = size(Q,1);
[UA,LambdaA]=eigs(A,r);
[UQ,LambdaQ]=eigs(Q,4);
%UA = sparse(UA);
%UQ = sparse(UQ);
H=H./sum(sum(H));y = H(:);

%tic;
% compute D 
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
t1 = toc;
end

tic;
D=1./sqrt(deg);
D(D==Inf)=0;
D = max(max(D))*ones(size(D, 1), 1);
X = D.*N;
% X=N;
iD2 = X.^2;
%iD2 = sparse(diag(iD2));
iX = 1./X;
iX(iX==Inf)=0;
%iX = sparse(diag(iX));


% X = sparse(diag(X));
% iX = 1./X;
% iX(iX==Inf)=0;
% iX=sparse(iX);
% I = sparse(diag(ones(m*n,1)));
% tic;
% D2 = (iX)^2;
% iD2 = sparse(1./D2); % can be trasformed to single
% iD2(iD2==Inf)=0;

clear D2;
%UQ=single(UQ);% can be trasformed to single
%UA=single(UA);% can be trasformed to single
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
%V2 = kron(UQ', UA'); %V2=U2';
%eta = U2'*bsxfun(@times, iD2, U2);

eta = U2'*bsxfun(@times,iD2,U2); %% time consuming!!!
temp = sparse(pinv(pinv(C2)-alpha.*eta));
clear eta;
clear C2;
%iD2=sparse(diag(iD2));
%Y = alpha.*(U2*(temp*(U2'*iD2')));%%
%Y = iD2*Y;
iXD2 = iX.*iD2;
s = (1-alpha).*(y+alpha.*(iXD2.*(U2*(temp*(U2'*(iXD2.*y))))));
clear temp;
clear iD2;
clear U2;
%Y = double(Y);
%s = (1-alpha).*(I + iX*Y*iX)*y;
t2=toc;
if issecondrun==1
    fprintf('running time = %f\n', t1+t2);
    time = t1+t2;
else
    fprintf('running time = %f\n', t2);
    time = t2;
end
SS=reshape(s,m,n);
end