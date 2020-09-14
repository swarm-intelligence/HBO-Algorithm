function [Out, nfe, fitness_result] = HBO_Sort(X,s,size,nfe,dim,max_nfe)
F_result = zeros(1, s);
for i = 1:s
    F_result(1,i) = F2(X(i,:));
    nfe = nfe + 1;
    if(nfe >= max_nfe)
       break;
    end
end

[~, sorted_index] = sort(F_result,2);

sorted = zeros(size,dim);

for j=1:size
    sorted(j,:) = X(sorted_index(1,j),:);
end
Out = sorted;
fitness_result = F_result;


