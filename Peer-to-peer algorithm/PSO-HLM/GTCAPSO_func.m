% function [gbest,gbestval,fitcount,suc,suc_fes]= MQ_PSO4_func(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax,varargin)
function [jg,gbest,gbestval,fitcount,suc,suc_fes,diversity,everyfit]= GTCAPSO_func(jingdu,fhd,Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax,varargin)
global orthm best_f best_keep initial_flag fbias gbias    
%   fprintf(' jingdu=');
%   jingdu
%   fprintf('\n Max_Gen=');
%   Max_Gen
%    fprintf('\n Particle_Number=');
%    Particle_Number
%    fprintf('\n Dimension=');
%    Dimension
%    fprintf('\n VRmin=');
%    VRmin
%    fprintf('\n VRmax=');
%    VRmax
%    fprintf('\n varargin=');
%    varargin
   
rand('state',sum(100*clock));
me=Max_Gen+1; 
ps=Particle_Number; 
original_ps=ps; 
D=Dimension;
iwt=0.7298*ones(1,D);
cc=[1.0, 1.0, 1.0];  
if length(VRmin)==1 
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
pos=VRmin+(VRmax-VRmin).*rand(ps,D);
vel=VRmin+2.*VRmax.*rand(ps,D);
mv=0.2.*(VRmax-VRmin); 

neighbor(1,:)=[ps,2]; 
for i=2:(ps-1)
    neighbor(i,:)=[i-1,i+1];
end
neighbor(ps,:)=[ps-1,1];
%  fprintf('22222222222222222222222222222222\n');
%  pos(1,:)
for i=1:ps  
	e(i,1)=feval(fhd,pos(i,:)',varargin{:}); 
    pbestval(i,1)=e(i,1);  
end
fitcount=ps; 
old_e=e; old_pos=pos;  
pbest=pos; 
old_pbest=pos;   
old_pbestval=pbestval;  


ratio=0.5;       
upper_elite=ps;      
lower_elite=round(ps/2);      

elites=round(ps/4);
sizeof_Af=round(ps/4);
sizeof_Adf=round(ps/4);
sizeof_Adfdd=round(ps/4);

archive_maxsize = ps; 

[gbestval,gbestid]=min(pbestval);   
gbest=pbest(gbestid,:);

for k=1:ps
    [~,tmpid]=min(pbestval(neighbor(k,:)));  
     aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(pbest(neighbor(k,tmpid),:)-pos(k,:)); 
    vel(k,:)=iwt(1,:).*vel(k,:)+aa(k,:);
    vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    pos(k,:)=(pos(k,:)>VRmax(k,:)).*VRmax(k,:)+(pos(k,:)<=VRmax(k,:)).*pos(k,:);
    pos(k,:)=(pos(k,:)<VRmin(k,:)).*VRmin(k,:)+(pos(k,:)>=VRmin(k,:)).*pos(k,:);
    e(k,1)=feval(fhd,pos(k,:)',varargin{:});
    fitcount=fitcount+1;       
    jg(1:fitcount)=gbestval;  
    tmp=(pbestval(k)<e(k));
    temp=repmat(tmp,1,D);  
    pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:); 
    pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);
    if pbestval(k)<gbestval
        gbest=pbest(k,:);
        gbestval=pbestval(k);
        gbestrep=repmat(gbest,ps,1); 
    end   
end 

tmppbest=[old_pbest;pbest]; 
tmppbestval=[old_pbestval; pbestval]; 
[tmppbestval, IX]=sort(tmppbestval); 
tmppbest=tmppbest(IX(:),:);         
Q_f(1)=tmppbestval(1); 
Q_f_pos(1,:)=tmppbest(1,:); 
index=1; 

for i=2:2*ps
    if ismember(tmppbest(i,:),Q_f_pos,'rows')==0  
        index=index+1; 
        Q_f(index,1)=tmppbestval(i);  
        Q_f_pos(index,:)=tmppbest(i,:);
        if index>=sizeof_Af            
           break; 
        end
    end     
end
prob_f=[index:-1:1];   
prob_f=prob_f/sum(prob_f);   

pbestval_dif=old_pbestval - pbestval;              
[pbestval_dif, IX]=sort(pbestval_dif,'descend'); 
tmppbest_dif=pbest(IX(:),:);      

Q_df(1)=pbestval_dif(1);
Q_df_pos(1,:)=tmppbest_dif(1,:); 
index=1;  
for i=2:ps           
    if ismember(tmppbest_dif(i,:),Q_df_pos,'rows')==0  
        index=index+1;  
        Q_df(index,1)=pbestval_dif(i);   
        Q_df_pos(index,:)=tmppbest_dif(i,:);
        if index>=sizeof_Adf        
           break; 
        end
    end     
end
prob_df=[index:-1:1];  
prob_df=prob_df/sum(prob_df);   

df=old_pbestval - pbestval;  
dd=sqrt(sum((pbest(:,:)-old_pbest(:,:)).*(pbest(:,:)-old_pbest(:,:)),2)); 

tmp_df_dd=df./dd; 
[tmp_df_dd, IX]=sort(tmp_df_dd,'descend'); 
tmppos_df_dd=pos(IX(:),:);  
Q_df_dd(1)=tmp_df_dd(1); Q_df_dd_pos(1,:)=tmppos_df_dd(1,:); 
index=1;      
for i=2:ps          
    if ismember(tmppos_df_dd(i,:),Q_df_dd_pos,'rows')==0  
        index=index+1;   
        Q_df_dd(index,1)=tmp_df_dd(i);     
        Q_df_dd_pos(index,:)=tmppos_df_dd(i,:);
        if index>=sizeof_Adfdd
           break;        
        end
    end     
end
prob_df_dd=[index:-1:1];  
prob_df_dd=prob_df_dd/sum(prob_df_dd);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----
% old2 = 1; new2 = fitcount;  %%%------记录评价次数的改变
% de_dd_size=size(Q_df_dd_pos,1);
% center = mean(Q_df_dd_pos);     
% div_i = 0;      %%%
% for di=1:de_dd_size    
%     div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
% end             %%%将所有维度的位置的值相加
% diversity(old2:new2) = div_i/(de_dd_size*(VRmax(1,1)-VRmin(1,1))) ; %%% 求 多样性   d（j） =（ X ij - 平均（X j））1到m的累加和
% %%%%%%%%%%%%%%%%%%

old_e=e; old_pos=pos;        %%% 将现在的位置变为上一代的位置，以便后续计算 位置变化
old_pbestval=pbestval; old_pbest = pbest;       %%% 将现在的适应度变为上一代的适应度，以便后续计算 适应度变化

recorded = 0;  % 达到精度时记录相关信息
suc = 0;        %%%-----
suc_fes = 0;    %%%-----

%%%% 以下为画图用
old = 1;
new = fitcount;
yyy(old:new) = gbestval;
old = new;
iter=1;
ceased=0;

archieve_val=[];%%%初始化优秀粒子 适应度 档案
archieve_pos=[];%%%初始化优秀粒子 位置 档案
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jg(fitcount)=gbestval;%最后结果保存
% while fitcount<Max_FES
while iter<me
    diversity(iter)=sum(sqrt(sum((pos-mean(pos)).^2,2)))/ps;
    everyfit(iter)=gbestval;
    iter=iter+1;  %%%  iter 先+1，因为在刚开始就已经更新过一次  所以从第二次开始 
%%%-------------- 培养优秀粒子-------------
    for k=1:ps
        % 选择学习目标个体          
        [obj_x(k,:),prob_f,prob_df,prob_df_dd]=select(iter,me,k,ps,sizeof_Af, Q_f,Q_f_pos,prob_f, Q_df,Q_df_pos,prob_df, Q_df_dd,Q_df_dd_pos,prob_df_dd, VRmax,VRmin);       
%        此时 obj_x 可进行评价，并根据结果将其纳入到种群中去
        obj_val(k)=feval(fhd,obj_x(k,:)',varargin{:});  %%%对其进行测试评价
        if iter>1/3*me
            obj_val1(k)=obj_val(k)*abs(randn);
        else
            obj_val1(k)=obj_val(k);
        end
        fitcount=fitcount+1;    %%%评价完之后 将其+1，并记录++++++++++++++++++++++++
         jg(fitcount)=gbestval;%结果保存
    end 
%%   
    [~, IX]=sort(obj_val1); 
       for k=1:ps
           j=IX(k); %%% 从IX（1）开始选择   即从 gbest 逐渐降低 开始选择
        if obj_val(j)<=gbestval  % 生成的榜样比 gbest 好
            [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x(j,:), obj_val(j));
            iwt=0.45 + 0.005 * tan(pi * (rand(1, D) - 0.5));        %%%计算保留系数
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(3).*rand(1,D).*(obj_x(j,:)-pos(j,:));   % （3）         %%% 按照 扩散型进行速度更新
            gbest(1,:)=obj_x(j,:);
            gbestval=obj_val(j);
            
         elseif obj_val1(j)<=gbestval %nbval% 生成的榜样比其邻居好  %生成的榜样比 pbest 好 
               iwt=0.45 + 0.005 * tan(pi * (rand(1, D) - 0.5));        %%%计算保留系数
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(3).*rand(1,D).*(obj_x(j,:)-pos(j,:))+0.001*abs(randn);   % （3）         %%% 按照 扩散型进行速度更新
            
        elseif obj_val1(j)<=pbestval(j)  %nbval% 生成的榜样比其邻居好  %生成的榜样比 pbest 好 
           [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x(j,:), obj_val(j));
            iwt=0.7298 + 0.005 * tan(pi * (rand(1, D) - 0.5));       %%%计算保留系数
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(obj_x(j,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  （2）%%% 按照 自信型进行速度更新
           
        elseif obj_val1(j)>pbestval(j)  % nbval % 生成的榜样比其邻居差  %生成的榜样比 pbest 差     
             iwt=0.7298+ 0.005 * tan(pi * (rand(1, D) - 0.5));         %%%计算保留系数
             vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(pbest(j,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  （1）%%% 按照 温和型进行速度更新    
        end
            
        vel(j,:)=(vel(j,:)>mv(1,:)).*mv(1,:)+(vel(j,:)<=mv(1,:)).*vel(j,:);     %%% 限制个体速度
        vel(j,:)=(vel(j,:)<(-mv(1,:))).*(-mv(1,:))+(vel(j,:)>=(-mv(1,:))).*vel(j,:);  %%%限制个体速度
        
        % 更新个体位置
        pos(j,:)= pos(j,:)+vel(j,:);  %%%更新个体位置
        pos(j,:)=(pos(j,:)>VRmax(1,:)).*(VRmax(1,:))+(pos(j,:)<=VRmax(1,:)).*pos(j,:);   %%% 限制个体位置
        pos(j,:)=(pos(j,:)<VRmin(1,:)).*(VRmin(1,:))+(pos(j,:)>=VRmin(1,:)).*pos(j,:);   %%%限制个体位置
        
         % 越界处理
        if rand<0.5
            pos(j,:)=(pos(j,:)>VRmax(1,:)).*VRmax(1,:)+(pos(j,:)<=VRmax(1,:)).*pos(j,:); 
            pos(j,:)=(pos(j,:)<VRmin(1,:)).*(VRmin(1,:))+(pos(j,:)>=VRmin(1,:)).*pos(j,:);
        else
            pos(j,:)=((pos(j,:)>=VRmin(1,:))&(pos(j,:)<=VRmax(1,:))).*pos(j,:)...
                +(pos(j,:)<VRmin(1,:)).*(VRmin(1,:)+0.1.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D))...
                +(pos(j,:)>VRmax(1,:)).*(VRmax(1,:)-0.1.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D));           
        end
        %评价个体
        e(j,1)=feval(fhd,pos(j,:)',varargin{:}); 
%         e(j,1)=feval(fhd,pos(j,:)',varargin{:})-fbias(func_num);
        fitcount=fitcount+1;        %%%+++++++++++++++++++++++++++++++++
         jg(fitcount)=gbestval;%结果保存
         
        if e(j,1)<pbestval(j,1)
           pbestval(j,1)=e(j,1);
           pbest(j,:)=pos(j,:);
           if pbestval(j,1)<gbestval
               gbestval=pbestval(j,1);
               gbest=pbest(j,:);
           end
        end
        
       end
       
      jg(fitcount)=gbestval;%最后结果保存
      
    %%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----        
%     old2 = new2; 
%     new2 = fitcount;
%     de_dd_size=size(Q_df_dd_pos,1);
%     center = mean(Q_df_dd_pos);
%     div_i = 0;
%     for di=1:de_dd_size    
%         div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
%     end
%     diversity(old2:new2) = div_i/(de_dd_size*(VRmax(1,1)-VRmin(1,1))) ; 
    %%%%%%%%%%%%%%%%%%

    %% Reusing Exemplars
%    % 将 archieve 中的个体与 pbest 比较和贪婪保存
  [pbest, pbestval, vel]=Reuse_Exemplars(pbest, pbestval, vel, archieve_pos, archieve_val);
   
   % 上述操作是将种群搜索历史过程中的个体历史最优解存入档案。
   % 有可能会将精英个体的多次个体历史最优解存入该档案，是否合理？如何解释
   % 此外，后面的 oldpbest 似乎也没有在每代中进行更新，需要修改部分代码
   % 先注释掉上面一句；然后在 Update_Q 函数后面添加 old_pbestval = pbestval 等几条语句
   
   
%%
    % 更新种群历史最优解
    [~, Index]=min(pbestval);
    if pbestval(Index)<gbestval 
        gbest=pbest(Index,:);
        gbestval=pbestval(Index);
%         gbestrep=repmat(gbest,ps,1);%update the gbest   
        ceased=0;   %%%-----------
    else
        ceased=ceased+1;        %%%----------------
    end  
   
      
    % 计算适应值变化信息 ------------   
     pbestval_df=old_pbestval - pbestval;  % 以 pbest 的变化为指标 
%      pbestval_df=old_e - e;  % 以 pos 的变化为指标
    
%     [pbestval_dif IX]=sort(pbestval_dif,'descend');  % 适应值降序排列
 
    % 计算 适应值变化/距离变化 信息
%     df=old_e-e; % 适应值变化
%     dd=sqrt(sum((pos(:,:)-old_pos(:,:)).*(pos(:,:)-old_pos(:,:)),2)); % 位置变化
    
    df=old_pbestval - pbestval;   % 使用 pbestval 
    dd=sqrt(sum((pbest(:,:)-old_pbest(:,:)).*(pbest(:,:)-old_pbest(:,:)),2)); % pbest位置变化
    
    df_dd=df./dd; % 搜索效率
%     [tmp_df_dd IX]=sort(tmp_df_dd,'descend');  % 搜索效率降序排列
%     tmppos_df_dd=pos(IX(:),:);          % 对应个体也排好序
  
    
    % 更新 Q_f,Q_df,Q_df_dd 及每个队列中个体的优先权
    [Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd]=...
        Update_Q(Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd, pbestval,pbest,  pbestval_df, df_dd,pos);
      
    % 新添加的
    old_pbestval = pbestval; old_pbest = pbest;
%     old_e = e; old_pos = pos;
    
%     end

%     if mod(i,50)==0,
%     plot(pos(:,D-1),pos(:,D),'b*');
%     hold on
%     plot(gbest(D-1),gbest(D),'r*');   
%     hold off
%     axis([VRmin(1,D-1),VRmax(1,D-1),VRmin(1,D),VRmax(1,D)])
%     title(['PSO: ',num2str(i),' generations, Gbestval=',num2str(gbestval)]);  
%     drawnow
% end 
%     if fitcount>=Max_FES
%         break;
%     end

end
%%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----        
%     old2 = new2; 
%     new2 = fitcount;
%     de_dd_size=size(Q_df_dd_pos,1);
%     center = mean(Q_df_dd_pos);
%     div_i = 0;
%     for di=1:de_dd_size    
%         div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
%     end
%     diversity(old2:new2) = div_i/(de_dd_size*(VRmax(1,1)-VRmin(1,1))) ; 
    %%%%%%%%%%%%%%%%%%
    
% diversity(Max_FES)=diversity(new2);
% xlabel('fes');
% ylabel('diversity');
% set(gca, 'Fontname', 'Times New Roman','FontSize',9);
% x=1:round(Max_FES/100):Max_FES;
% x(101)=Max_FES;
% diversity(Max_FES) = gbestval;
% hold on;
% plot(x,diversity(x),'-r','MarkerFaceColor','r','MarkerSize',5);  

end

