
%%%%%%%%%%%%%%%%%% �������õ������¸���� pbest  %%%%%%%%%%%%%%%%%%%%%%
function [pbest, pbestval, vel]=Reuse_Exemplars(pbest, pbestval, vel, archieve_pos, archieve_val)
    archiev_size=length(archieve_val); %%% ������Ao����������
    [ps,D]=size(pbest); %%% pbest��С�����ڱ���
   
    [Val_p, Index_p]=sort(pbestval,'descend');%%%��pbest��Ӧ������ ��������Val_p��
    i_p=1;
    
    index=1;
    rp=randperm(archiev_size);%%%����������ҷֲ�
    for i=1:floor(archiev_size/2) %%% ���з���---2��Ϊһ�飬Ѱ�����Ž��б��浽tmppos��tmpval��
        if archieve_val(rp(2*i-1))<archieve_val(rp(2*i))%%% ���з���---2��Ϊһ�飬Ѱ�����Ž��б��浽tmppos��tmpval��
           tmpval(index)=archieve_val(rp(2*i-1));
           tmppos(index,:)=archieve_pos(rp(2*i-1),:);
           index=index+1;
        else
           tmpval(index)=archieve_val(rp(2*i));%%% ���з���---2��Ϊһ�飬Ѱ�����Ž��б��浽tmppos��tmpval��
           tmppos(index,:)=archieve_pos(rp(2*i),:); 
           index=index+1;
        end
    end
    
    maxsize=length(tmpval);
    index=1;
    while i_p<=ps && index<=maxsize   %%%���б����滻
        if Val_p(i_p,1)>tmpval(index) && rand<0.5
           pbestval(Index_p(i_p))=tmpval(index);
           pbest(Index_p(i_p),:)=tmppos(index,:);  

           vel(Index_p(i_p),:)=zeros(1,D);    
           i_p=i_p+1;
        end        
        index=index+1;
    end

end
