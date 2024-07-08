
%%%%%%%%%%%%%%%%%% 重新利用档案更新个体的 pbest  %%%%%%%%%%%%%%%%%%%%%%
function [pbest, pbestval, vel]=Reuse_Exemplars(pbest, pbestval, vel, archieve_pos, archieve_val)
    archiev_size=length(archieve_val); %%% 论文中Ao档案的容量
    [ps,D]=size(pbest); %%% pbest大小，用于遍历
   
    [Val_p, Index_p]=sort(pbestval,'descend');%%%将pbest适应度排序 并保存在Val_p中
    i_p=1;
    
    index=1;
    rp=randperm(archiev_size);%%%将其随机打乱分布
    for i=1:floor(archiev_size/2) %%% 进行分组---2个为一组，寻找最优进行保存到tmppos和tmpval中
        if archieve_val(rp(2*i-1))<archieve_val(rp(2*i))%%% 进行分组---2个为一组，寻找最优进行保存到tmppos和tmpval中
           tmpval(index)=archieve_val(rp(2*i-1));
           tmppos(index,:)=archieve_pos(rp(2*i-1),:);
           index=index+1;
        else
           tmpval(index)=archieve_val(rp(2*i));%%% 进行分组---2个为一组，寻找最优进行保存到tmppos和tmpval中
           tmppos(index,:)=archieve_pos(rp(2*i),:); 
           index=index+1;
        end
    end
    
    maxsize=length(tmpval);
    index=1;
    while i_p<=ps && index<=maxsize   %%%进行遍历替换
        if Val_p(i_p,1)>tmpval(index) && rand<0.5
           pbestval(Index_p(i_p))=tmpval(index);
           pbest(Index_p(i_p),:)=tmppos(index,:);  

           vel(Index_p(i_p),:)=zeros(1,D);    
           i_p=i_p+1;
        end        
        index=index+1;
    end

end
