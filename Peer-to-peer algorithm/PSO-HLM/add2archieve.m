
%%%%%%%%%%%%%%%%%% ���µ���  %%%%%%%%%%%%%%%%%%%%%%
function [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, max_size, obj_x, obj_val)
    archiev_size=length(archieve_val); %%%�����У��������ӣ�Ao����������
    if archiev_size<max_size  % ��������δ�ﵽ����ֵ��ֱ�ӷ��� %
        archiev_size=archiev_size+1;
        archieve_val(archiev_size,1)=obj_val;
        archieve_pos(archiev_size,:)=obj_x;   
    else  % ���������Ѵﵽ����ֵ���滻���ĸ��� %
%         [tmpB, tmpIX]=max(archieve_val);
%         archiev_size=length(archieve_val);
% %         if tmpB>obj_val  % obj_x�ȵ�ǰ������������Ҫ��ʱ���滻�ø���
%         archieve_val(tmpIX,1)=obj_val;
%         archieve_pos(tmpIX,:)=obj_x;
% %         end
        % ���������Ѵﵽ����ֵ���ڵ�������ѡ2�����滻���нϲ�ĸ��� %  %%% ������
        rp=randperm(archiev_size);
        if archieve_val(rp(1))<archieve_val(rp(2))% && archieve_val(rp(2))>obj_val % �����ˣ�ע�⣡����������        
            archieve_val(rp(2),1)=obj_val;
            archieve_pos(rp(2),:)=obj_x;
                        
        else%if archieve_val(rp(2))<archieve_val(rp(1)) && archieve_val(rp(1))>obj_val %  % �����ˣ�ע�⣡����������
                archieve_val(rp(1),1)=obj_val;
                archieve_pos(rp(1),:)=obj_x;            
        end

    end
end