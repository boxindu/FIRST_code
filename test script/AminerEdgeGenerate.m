load('aminer_chenchen.mat');
n = size(author_author,1);
x = size(paper_authors,1);
edge_label2=sparse(1:n,1:n,0,n,n);
perm_result = [];
for q = 1:x
    index = find(paper_authors(q,:));
    l=size(index,2);
    index = cat(1, index, index);
    total = size(index,2)^2;
    for i = 0:total - 1
    val = i;
    info = zeros(1,2);
    for j = 2:-1:1
        info(j) = mod(val,l);
        val = fix(val./l);
    end
    tem = [];
    for k = 0:2-1
        tem = [tem, index(k+1, info(k+1)+1)];
    end
    perm_result = [perm_result; tem];
    end
end
for i = 1:size(perm_result,1)
    in = perm_result(i,:);
    edge_label2(in(1),in(2)) = edge_label2(in(1),in(2))+1;
end


    