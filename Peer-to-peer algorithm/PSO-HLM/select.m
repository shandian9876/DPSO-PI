
%%%%%%%%%%%%%%%%%% 为子种群每个个体选择合适的学习对象  %%%%%%%%%%%%%%%%%%%%%%
function [obj_x,prob_f,prob_df,prob_df_dd]=select(iter,me,k,ps,elites,Queue_f,pos_f,prob_f, Queue_df,pos_df,prob_df, Queue_df_dd, pos_df_dd, prob_df_dd,  Xmax,Xmin)

    [obj_x, prob_f, prob_df]=Generate(iter,me,Queue_f,prob_f,pos_f,  Queue_df_dd,prob_df_dd,pos_df_dd, Xmax,Xmin);  % 

end
