%% 计算Novelty函数
%==========================================================================
function Novelty = Computing_Novelty(NP, D, pbest, gbest, num_neighbors, type)
    switch type
    case 1
    % 个体的 pbest 与最临近的 num_neighbors 个个体 pbest 的距离来表征 Novelty
         for i = 1:NP
            for j = 1:NP
                % %欧几里得距离 
%                distance(1,j) = norm(pbest(i,:)-pbest(j,:));
                % %曼哈顿距离 
%                 distance(1,j) = pdist([pbest(i,:); pbest(j,:)], 'cityblock');%'cityblock'表示计算曼哈顿距离。 
                  distance(1,j) = sum(abs(pbest(i,:)-pbest(j,:)));
            end        

            distance(distance==0)=[]; % 删掉 0 元素，即不考虑个体与其自身的距离
            if ~isempty(distance)
                [Novel_val, Novel_IX] = sort(distance,'ascend');

                if length(Novel_IX) == 1   % 当种群已经收敛时:对 pbest 随机扰动
%                     Novelty(i,1) = mean(Novel_val(Novel_IX(1)));
%                     sigma = (lu(2)-lu(1))/200;
%                     pbest(i,:) = pbest(i,:) + normrnd(0,sigma,1,D);

                else
                    Novelty(i,1) = mean(Novel_val(Novel_IX(1:min(num_neighbors,length(Novel_IX)))));
                end               
            else
                Novelty(i,1) = 0;
            end
         end
     
         if max(Novelty) > 0
            Novelty = Novelty/max(Novelty);     % 归一化，值越大表明越新颖！
         else
            return;
         end
    end
end