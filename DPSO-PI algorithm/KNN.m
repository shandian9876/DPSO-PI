function [index_number] = KNN(pbest)
% function [ODA_AbnormalObject_Number,ODA_NormalObject_Number] = KNN(pbest)

x=pbest;


[m,n]=size(x);
k=10;%近邻个数
Abnormal_number=2;%离群对象个数

Dist = calculate_manhattan_distances(x);
SortDist=sort(Dist,2,'ascend');%2代表数据按行进行排序
Nei_k=SortDist(:,1:k+1);%因为对象离其自身的距离为0，所以再多考虑一个对象
OF=1./sum(Nei_k,2);

[OF_value,index_number]=sort(OF,'descend');
% ODA_AbnormalObject_Number=(index_number(m-Abnormal_number+1:end,:))';%outlier detection algorithm 算法认定的异常对象的编号
% ODA_NormalObject_Number=index_number(1:2,:)';%outlier detection algorithm算法认定的正常对象的编号

end


function distances = calculate_manhattan_distances(pbest)
    [num_points, dim] = size(pbest);
    distances = zeros(num_points, num_points);
    
    for i = 1:num_points
        for j = i+1:num_points
            distances(i, j) = sum(abs(pbest(i, :) - pbest(j, :)));
            distances(j, i) = distances(i, j); % Since it's symmetric, we can assign the same value to distances(j, i)
        end
    end
end

