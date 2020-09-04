%% This demo is using Gordian_v2 dataset

clear all;

addpath('./data');
load('C:./data/Gordian_V2.mat');
S = FINAL(graph, signal, [], [], {}, {}, ones(1003,19)/(1003*19), 0.5, 40, 1e-10);
[~,Sindex] = sort(S,2,'descend');
temp2 = Sindex(:,1:20);
temp2 = unique(temp2);
signal = signal +1;
num=0;
for i = 1:19
if ismember(signal(i), temp2)
num=num+1;
end
end

Precision = num/length(temp2);
Recall = num/length(signal);
rmpath('./data');

%% visualization

Gord=digraph(graph);
Sig = subgraph(Gord, signal);
figure('Name', 'Signal'); plot(Sig);
First=subgraph(Gord, temp2);
figure('Name', 'matching subgraph'); plot(First);

