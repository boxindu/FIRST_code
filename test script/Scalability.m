% load('DBLP.mat');
TIME = [];
TIME2=[];
DIFF=[];
DIST=[];
label = node_label2;
for i = 1:20
    j =i;
    Q1=Q;
%    Ed = edge_label2;
    while j>0
        Q1 = cat(2, Q1,Q);
%        edge_label2=cat(2,edge_label2,Ed);
        j=j-1;
    end
    j = i;
    Q2=Q1;
%    Ed = edge_label2;
    while j>0
        Q1 = cat(1, Q1,Q2);
%        edge_label2 = cat(1, edge_label2, Ed);
        j=j-1;
    end
    
node_label2=cat(1,node_label2,label);
H=ones(13*(i+1),9143)./(13*(i+1)*9143);
% H=ones(5*(i+1),9143)./(5*(i+1)*9143);
max_iter=50;
[~, time1] = final_N(A, Q1, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
[~, time2] = finalN_low(A, Q1, H, node_label1, node_label2, l_node, alpha, r, relax);
[~, time3] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, r, relax, 1);
[~, time4] = INCREMENTAL_N(A, Q1, H, node_label1, node_label2, node_label2, l_node, alpha,r, relax, 1);
% [~,time5]=final_NE(A, Q1, H, node_label1, node_label2, edge_label1, edge_label2, l_node, l_edge, alpha, max_iter, relax);
% [~,time6] = INCREMENTAL_E(A, Q1, H, node_label1, node_label2, l_node, edge_label1, ...
% edge_label2, l_edge, alpha, r, relax, 1, 1,2, 0);
TIME = [TIME; time1, time2, time3, time4];
% TIME2 = [TIME2;time5,time6];

% [~, time1] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, 2, relax, 1);
% [~, time2] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha,3, relax, 1);
% [~, time3] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, 5, relax, 1);
% [~, time4] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, 8, relax, 1);
% [~, time5] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, 10, relax, 1);
% TIME2=[TIME2;time1,time2,time3,time4,time5];

% [SS,~] = INCREMENTAL_Q(A, Q1, H, node_label1, node_label2, l_node, alpha, r, relax, 1);
% [SN, ~] = INCREMENTAL_N(A, Q1, H, node_label1, node_label2, node_label2, l_node, alpha,r, relax, 1);
% [S, time] = final_N(A, Q1, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
% diff1=abs(norm(SS,'fro')-norm(S,'fro'));
% dist1 = diff1/abs(norm(S,'fro'));
% diff2=abs(norm(SN,'fro')-norm(S,'fro'));
% dist2 = diff2/abs(norm(S,'fro'));
% DIFF = [DIFF;diff1,diff2];
% DIST = [DIST; dist1,dist2];
end