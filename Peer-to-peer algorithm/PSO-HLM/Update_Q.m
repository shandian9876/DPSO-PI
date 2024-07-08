function [Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd]=Update_Q(Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd, ...
                 pbestval,pbest,  pbestval_df, df_dd,pos)
    % ԭ������Ϣ
%     [B_f IX_f]=max(Q_f);
%     [B_df IX_df]=min(Q_df);
%     [B_df_dd IX_df_dd]=min(Q_df_dd);
%     
%     [bestv  bestIX]=sort(pbestval);    
    
    %--------- �����ǽ�������Ⱥ���ڸ�ָ�������ŵķ�����Ӧ���� -------------%
    %--------- �Ƿ�Ӧ�ý�������Ⱥ�и�ָ�������ŵķ�����Ӧ���� ������-------%
    
    % ��Ӧֵ���� ----------------------   
    [tmppbestval, IX]=sort(pbestval);  % ��Ӧֵ����������
    tmppbest=pbest(IX(:),:);          % ��Ӧ����Ҳ�ź���
    ps=length(Q_f);
%     tmpQ_f=[];
%     tmpQ_f_pos=[];
    tmpQ_f=[tmppbestval(1)];
    tmpQ_f_pos=[tmppbest(1,:)];
    index=1; i=1; j=1;
    while index<=ps && i<=ps && j<=ps
        if tmppbestval(i)<Q_f(j,1)
            if ismember(tmppbest(i,:),tmpQ_f_pos,'rows')==0  % �� i ��Ԫ�ز��ڶ�����                
                tmpQ_f(index,1)=tmppbestval(i); 
                tmpQ_f_pos(index,:)=tmppbest(i,:); 
                index=index+1;
            end   
            i=i+1;             
        else 
            if ismember(Q_f_pos(j,:),tmpQ_f_pos,'rows')==0  % �� i ��Ԫ�ز��ڶ�����                
                tmpQ_f(index,1)=Q_f(j,1); 
                tmpQ_f_pos(index,:)=Q_f_pos(j,:); 
                index=index+1;                
            end 
            j=j+1;
        end
    end
    
%     Q_f=tmpQ_f;
%     Q_f_pos=tmpQ_f_pos;
    tmpsize=length(tmpQ_f);
    for i=1:tmpsize
        Q_f(i,1)=tmpQ_f(i,1);
        Q_f_pos(i,:)=tmpQ_f_pos(i,:); 
    end
    
    % ��Ӧֵ�仯 ���� ------------------------------
    [tmppbestval_df, IX]=sort(pbestval_df,'descend');  % ��Ӧֵ����������
    tmppbest_df=pbest(IX(:),:);          % ��Ӧ����Ҳ�ź���
    ps=length(Q_df);
    Q_df=tmppbestval_df(1:ps,:);
    Q_df_pos=tmppbest_df(1:ps,:); 
%     for i=1:ps
%         Q_df(i,1)=tmppbestval_df(i,1);
%         Q_df_pos(i,:)=tmppbest_df(i,:); 
%     end
    
    % ��Ӧֵ�仯/���о���仯 ���� ------------------
    [tmp_df_dd, IX]=sort(df_dd,'descend');  % ��Ӧֵ����������
    tmp_pos=pos(IX(:),:);          % ��Ӧ����Ҳ�ź���
    ps=length(Q_df_dd);
%     tmpQ_df_dd=[];
%     tmpQ_df_dd_pos=[];
%     index=1; i=1; j=1;
%     while index<=ps && i<=ps && j<=ps
%         if tmpdf_dd(i)>Q_df_dd(j,1)
%             if ismember(tmpdf_dd(i),tmpQ_df_dd)==0  % �� i ��Ԫ�ز��ڶ�����                
%                 tmpQ_df_dd(index,1)=tmpdf_dd(i); 
%                 tmpQ_df_dd_pos(index,:)=tmppos(i,:); 
%                 index=index+1;
%             end   
%             i=i+1;             
%         else 
%             if ismember(Q_df_dd(j,1),tmpQ_df_dd)==0  % �� i ��Ԫ�ز��ڶ�����                
%                 tmpQ_df_dd(index,1)=Q_df_dd(j,1); 
%                 tmpQ_df_dd_pos(index,:)=Q_df_dd_pos(j,:); 
%                 index=index+1;                
%             end 
%             j=j+1;
%         end
%     end

    Q_df_dd=tmp_df_dd(1:ps,:);
    Q_df_dd_pos=tmp_pos(1:ps,:); 
%     tmpsize=length(tmp_df_dd);
%     for i=1:tmpsize
%         Q_df_dd(i,1)=tmp_df_dd(i,1);
%         Q_df_dd_pos(i,:)=tmp_pos(i,:); 
%     end
    
    
end    
