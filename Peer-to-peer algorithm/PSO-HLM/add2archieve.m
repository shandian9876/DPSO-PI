
%%%%%%%%%%%%%%%%%% 更新档案  %%%%%%%%%%%%%%%%%%%%%%
function [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, max_size, obj_x, obj_val)
    archiev_size=length(archieve_val); %%%论文中（优秀粒子）Ao档案的容量
    if archiev_size<max_size  % 档案容量未达到上限值，直接放入 %
        archiev_size=archiev_size+1;
        archieve_val(archiev_size,1)=obj_val;
        archieve_pos(archiev_size,:)=obj_x;   
    else  % 档案容量已达到上限值，替换最差的个体 %
%         [tmpB, tmpIX]=max(archieve_val);
%         archiev_size=length(archieve_val);
% %         if tmpB>obj_val  % obj_x比当前档案中最差个体要优时，替换该个体
%         archieve_val(tmpIX,1)=obj_val;
%         archieve_pos(tmpIX,:)=obj_x;
% %         end
        % 档案容量已达到上限值，在档案中任选2个，替换其中较差的个体 %  %%% 更改了
        rp=randperm(archiev_size);
        if archieve_val(rp(1))<archieve_val(rp(2))% && archieve_val(rp(2))>obj_val % 更改了，注意！！！！！！        
            archieve_val(rp(2),1)=obj_val;
            archieve_pos(rp(2),:)=obj_x;
                        
        else%if archieve_val(rp(2))<archieve_val(rp(1)) && archieve_val(rp(1))>obj_val %  % 更改了，注意！！！！！！
                archieve_val(rp(1),1)=obj_val;
                archieve_pos(rp(1),:)=obj_x;            
        end

    end
end