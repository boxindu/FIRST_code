% load('DBLP.mat');
n=13;
P=[];
C=[];
for k=1:13
    [SS, ~] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,2, relax, 1);
     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
%     [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
    [precision1, recall1] = AccuracyComp(S,SS,n,k);
    [SS, ~] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,3, relax, 1);
     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 3, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
    [precision2, recall2] = AccuracyComp(S,SS,n,k);
    [SS, ~] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,5, relax, 1);
     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 5, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
    [precision3, recall3] = AccuracyComp(S,SS,n,k);
    [SS, ~] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,8, relax, 1);
     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 8, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
    [precision4, recall4] = AccuracyComp(S,SS,n,k);
    [SS, ~] = INCREMENTAL_Q(A, Q, H, node_label1, node_label2, l_node, alpha,10, relax, 1);
     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 10, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
    [precision5, recall5] = AccuracyComp(S,SS,n,k);
    P = [P; precision1, precision2, precision3,precision4,precision5];
    C = [C; recall1,recall2,recall3,recall4,recall5];
end


% n=13;
% P=[];
% C=[];
% for k=1:13
%     [SS, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
% %     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
%     [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
%     [precision1, recall1] = AccuracyComp(S,SS,n,k);
%     [SS, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
% %     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 3, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
%     [precision2, recall2] = AccuracyComp(S,SS,n,k);
%     [SS, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
% %     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 5, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
%     [precision3, recall3] = AccuracyComp(S,SS,n,k);
%     [SS, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
% %     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 8, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
%     [precision4, recall4] = AccuracyComp(S,SS,n,k);
%     [SS, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 2, relax);
% %     [S, ~] = finalN_low(A, Q, H, node_label1, node_label2, l_node, alpha, 10, relax);
% [S, ~] = final_N(A, Q, H, node_label1, node_label2, l_node, alpha, max_iter, relax);
%     [precision5, recall5] = AccuracyComp(S,SS,n,k);
%     P = [P; precision1, precision2, precision3,precision4,precision5];
%     C = [C; recall1,recall2,recall3,recall4,recall5];
% end
