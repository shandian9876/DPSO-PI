function [Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd]=Update_Q(Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd, ...
                 pbestval,pbest,  pbestval_df, df_dd,pos)
    % 原队列信息
%     [B_f IX_f]=max(Q_f);
%     [B_df IX_df]=min(Q_df);
%     [B_df_dd IX_df_dd]=min(Q_df_dd);
%     
%     [bestv  bestIX]=sort(pbestval);    
    
    %--------- 这里是将整改种群中在各指标上最优的放入相应队列 -------------%
    %--------- 是否应该将各子种群中个指标上最优的放入相应队列 ？？？-------%
    
    % 适应值队列 ----------------------   
    [tmppbestval, IX]=sort(pbestval);  % 适应值升序序排列
    tmppbest=pbest(IX(:),:);          % 对应个体也排好序
    ps=length(Q_f);
%     tmpQ_f=[];
%     tmpQ_f_pos=[];
    tmpQ_f=[tmppbestval(1)];
    tmpQ_f_pos=[tmppbest(1,:)];
    index=1; i=1; j=1;
    while index<=ps && i<=ps && j<=ps
        if tmppbestval(i)<Q_f(j,1)
            if ismember(tmppbest(i,:),tmpQ_f_pos,'rows')==0  % 第 i 个元素不在队列中                
                tmpQ_f(index,1)=tmppbestval(i); 
                tmpQ_f_pos(index,:)=tmppbest(i,:); 
                index=index+1;
            end   
            i=i+1;             
        else 
            if ismember(Q_f_pos(j,:),tmpQ_f_pos,'rows')==0  % 第 i 个元素不在队列中                
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
    
    % 适应值变化 队列 ------------------------------
    [tmppbestval_df, IX]=sort(pbestval_df,'descend');  % 适应值升序序排列
    tmppbest_df=pbest(IX(:),:);          % 对应个体也排好序
    ps=length(Q_df);
    Q_df=tmppbestval_df(1:ps,:);
    Q_df_pos=tmppbest_df(1:ps,:); 
%     for i=1:ps
%         Q_df(i,1)=tmppbestval_df(i,1);
%         Q_df_pos(i,:)=tmppbest_df(i,:); 
%     end
    
    % 适应值变化/飞行距离变化 队列 ------------------
    [tmp_df_dd, IX]=sort(df_dd,'descend');  % 适应值升序序排列
    tmp_pos=pos(IX(:),:);          % 对应个体也排好序
    ps=length(Q_df_dd);
%     tmpQ_df_dd=[];
%     tmpQ_df_dd_pos=[];
%     index=1; i=1; j=1;
%     while index<=ps && i<=ps && j<=ps
%         if tmpdf_dd(i)>Q_df_dd(j,1)
%             if ismember(tmpdf_dd(i),tmpQ_df_dd)==0  % 第 i 个元素不在队列中                
%                 tmpQ_df_dd(index,1)=tmpdf_dd(i); 
%                 tmpQ_df_dd_pos(index,:)=tmppos(i,:); 
%                 index=index+1;
%             end   
%             i=i+1;             
%         else 
%             if ismember(Q_df_dd(j,1),tmpQ_df_dd)==0  % 第 i 个元素不在队列中                
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
