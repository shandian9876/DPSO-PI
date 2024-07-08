% function [gbest,gbestval,fitcount,suc,suc_fes]= MQ_PSO4_func(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax,varargin)
function [gbest,gbestval,fitcount,suc,suc_fes,diversity,everyfit]= TAPSO_func(jingdu,fhd,Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax,varargin)
    
    
global orthm best_f best_keep initial_flag fbias gbias                                                                
%[gbest,gbestval,fitcount]= PSO_local_func('f8',3500,200000,30,30,-5.12,5.12)
rand('state',sum(100*clock));
me=Max_Gen+1;
N=Particle_Number;
ps=Particle_Number;
original_ps=ps; 
D=Dimension;
everyfit=zeros(1,Max_Gen);

% cc=[2.05 2.05];   %acceleration constants
% iwt=0.9-(1:me).*(0.5./me);
% cc=[1.49445 1.49445];   %acceleration constants
iwt=0.6*ones(1,D);%0.7298;

% cc=[1.49445 1.49445];
 cc=[1.0, 1.0, 1.0];   %acceleration constants
%cc=[0.3, 0.6, 1.0];
% iwt=0.7;

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

% mv=0.2.*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
% Vmin=repmat(-mv,ps,1);
% Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);
vel=VRmin+2.*VRmax.*rand(ps,D);%initialize the velocity of the particles
mv=0.2.*(VRmax-VRmin);

neighbor(1,:)=[ps,2];
for i=2:(ps-1)
    neighbor(i,:)=[i-1,i+1];
end
neighbor(ps,:)=[ps-1,1];

for i=1:ps   % 记录每个分组的适应值 
	e(i,1)=feval(fhd,pos(i,:)',varargin{:});       
%    e(i,1)=feval(fhd,pos(i,:)',varargin{:})-fbias(func_num);
    pbestval(i,1)=e(i,1);    
end
fitcount=ps;
old_e=e; old_pos=pos;
pbest=pos; old_pbest=pos;
old_pbestval=pbestval;

% 加入队列的优势个体的比例
ratio=0.5;        % 正常情况下优势个体占整个种群的比例
upper_elite=ps;      % 优势个体数量的上限
lower_elite=round(ps/2);       % 优势个体数量的下限

elites=round(ps/4);  %round()取最近整数
sizeof_Af=round(ps/4);%
sizeof_Adf=round(ps/4);
sizeof_Adfdd=round(ps/4);
% elites=min(max(ceil(ps*ratio),lower_elite),upper_elite);

archive_maxsize = ps; % 用于记录生成的较优的榜样的档案大小

[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
% gbestrep=repmat(gbest,ps,1);

for k=1:ps
    [~,tmpid]=min(pbestval(neighbor(k,:)));
     aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(pbest(neighbor(k,tmpid),:)-pos(k,:));
%     vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); 
    vel(k,:)=iwt(1,:).*vel(k,:)+aa(k,:);
    vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    pos(k,:)=(pos(k,:)>VRmax(k,:)).*VRmax(k,:)+(pos(k,:)<=VRmax(k,:)).*pos(k,:); 
    pos(k,:)=(pos(k,:)<VRmin(k,:)).*VRmin(k,:)+(pos(k,:)>=VRmin(k,:)).*pos(k,:);

    e(k,1)=feval(fhd,pos(k,:)',varargin{:});
%    e(k,1)=feval(fhd,pos(k,:)',varargin{:})-fbias(func_num);
    fitcount=fitcount+1;
    tmp=(pbestval(k)<e(k));
    temp=repmat(tmp,1,D);
    pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
    pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
    if pbestval(k)<gbestval
        gbest=pbest(k,:);
        gbestval=pbestval(k);
        gbestrep=repmat(gbest,ps,1);%update the gbest  
    end   
end

%% ----------- 初始化 适应值 队列 -------------%
tmppbest=[old_pbest;pbest]; tmppbestval=[old_pbestval; pbestval];  % 以 pbest 的变化为指标
[tmppbestval, IX]=sort(tmppbestval);  % 适应值升序排列
% tmppbest=[old_pos;pos]; tmppbestval=[old_e; e];                      % 以 pos 的变化为指标
% [tmppbestval, IX]=sort(tmppbestval);  % 适应值升序排列

tmppbest=tmppbest(IX(:),:);          % 对应个体也排好序
Q_f(1)=tmppbestval(1); Q_f_pos(1,:)=tmppbest(1,:); % 将适应值及对应个体放入初始队列
index=1;

for i=2:2*ps
    if ismember(tmppbest(i,:),Q_f_pos,'rows')==0  % 第 i 个元素不在队列中
        index=index+1;
        Q_f(index,1)=tmppbestval(i); 
        Q_f_pos(index,:)=tmppbest(i,:);
        if index>=sizeof_Af              
           break; 
        end
    end     
end
prob_f=[index:-1:1];   
prob_f=prob_f/sum(prob_f);    % 队列中元素的候选概率

%% ----------- 初始化 适应值变化 队列 -----------%
% pbestval_dif=old_e - e;                                              % 以 pos 的变化为指标
% [pbestval_dif, IX]=sort(pbestval_dif,'descend');  % 适应值降序排列
% tmppbest_dif=pos(IX(:),:);          % 对应个体也排好序

pbestval_dif=old_pbestval - pbestval;                              % 以 pbest 的变化为指标
[pbestval_dif, IX]=sort(pbestval_dif,'descend');  % 适应值降序排列
tmppbest_dif=pbest(IX(:),:);          % 对应个体也排好序

Q_df(1)=pbestval_dif(1); Q_df_pos(1,:)=tmppbest_dif(1,:); % 将适应值及对应个体放入初始队列
index=1;
% elites=ceil(ps*ratio);
for i=2:ps
    if ismember(tmppbest_dif(i,:),Q_df_pos,'rows')==0  % 第 i 个元素不在队列中
        index=index+1;
        Q_df(index,1)=pbestval_dif(i); 
        Q_df_pos(index,:)=tmppbest_dif(i,:);
        if index>=sizeof_Adf
           break; 
        end
    end     
end
prob_df=[index:-1:1];   
prob_df=prob_df/sum(prob_df);    % 队列中元素的候选概率

%% ----------- 初始化 适应值变化/距离变化 队列 -----------%
% df=old_e - e;   % 使用 pbestval   
% dd=sqrt(sum((pos(:,:)-old_pos(:,:)).*(pos(:,:)-old_pos(:,:)),2)); % pos位置变化

df=old_pbestval - pbestval;   % 使用 pbestval 
dd=sqrt(sum((pbest(:,:)-old_pbest(:,:)).*(pbest(:,:)-old_pbest(:,:)),2)); % pbest位置变化

tmp_df_dd=df./dd; % 搜索效率
[tmp_df_dd, IX]=sort(tmp_df_dd,'descend');  % 搜索效率降序排列
tmppos_df_dd=pos(IX(:),:);          % 对应个体也排好序

Q_df_dd(1)=tmp_df_dd(1); Q_df_dd_pos(1,:)=tmppos_df_dd(1,:); % 将适应值及对应个体放入初始队列
index=1;
% elites=ceil(ps*ratio);
for i=2:ps
%     if ismember(tmp_df_dd(i),Q_df_dd,'rows')==0  % 第 i 个元素不在队列中
    if ismember(tmppos_df_dd(i,:),Q_df_dd_pos,'rows')==0  % 第 i 个元素不在队列中
        index=index+1;
        Q_df_dd(index,1)=tmp_df_dd(i); 
        Q_df_dd_pos(index,:)=tmppos_df_dd(i,:);
        if index>=sizeof_Adfdd
           break; 
        end
    end     
end
prob_df_dd=[index:-1:1];   
% prob_df_dd=prob_df_dd.*prob_df_dd;
prob_df_dd=prob_df_dd/sum(prob_df_dd);    % 队列中元素的候选概率

%%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----
old2 = 1; new2 = fitcount;
de_dd_size=size(Q_df_dd_pos,1);
center = mean(Q_df_dd_pos);
div_i = 0;
for di=1:de_dd_size    
    div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
end
% diversity(old2:new2) = div_i/(de_dd_size*(VRmax(1,1)-VRmin(1,1))) ; 
%%%%%%%%%%%%%%%%%%

old_e=e; old_pos=pos;
old_pbestval=pbestval; old_pbest = pbest;

recorded = 0;  % 达到精度时记录相关信息
suc = 0;
suc_fes = 0;

%%%% 以下为画图用
old = 1;
new = fitcount;
yyy(old:new) = gbestval;
old = new;
iter=1;
ceased=0;

archieve_val=[];
archieve_pos=[];

% while fitcount<Max_FES
 while iter<me
    
    diversity(iter)=sum(sqrt(sum((pos-mean(pos)).^2,2)))/N;
    everyfit(iter)=gbestval;
    iter=iter+1;    
    [~, IX]=sort(pbestval);

    for k=1:ps
        j=IX(k);
        % 选择学习目标个体          
        [obj_x,prob_f,prob_df,prob_df_dd]=select(k,ps,sizeof_Af, Q_f,Q_f_pos,prob_f, Q_df,Q_df_pos,prob_df, Q_df_dd,Q_df_dd_pos,prob_df_dd, VRmax,VRmin);
        
%        此时 obj_x 可进行评价，并根据结果将其纳入到种群中去
        obj_val=feval(fhd,obj_x(1,:)',varargin{:});
%         tmpval=feval(fhd,obj_x(1,:)',varargin{:})-fbias(func_num);
        fitcount=fitcount+1;

%%   
%         iwt=0.7298 + 0.005 * tan(pi * (rand(1, D) - 0.5));
        
        if obj_val<=gbestval  % 生成的榜样比 gbest 好
            [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x, obj_val);
            iwt=0.45 + 0.005 * tan(pi * (rand(1, D) - 0.5));
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(3).*rand(1,D).*(obj_x(1,:)-pos(j,:));   % （3）          
            gbest(1,:)=obj_x;
            gbestval=obj_val;
            
        elseif obj_val<=pbestval(j)  %nbval% 生成的榜样比其邻居好  %生成的榜样比 pbest 好 
            [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x, obj_val);
            iwt=0.6 + 0.005 * tan(pi * (rand(1, D) - 0.5));
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(obj_x(1,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  （2）
           
        elseif obj_val>pbestval(j)  % nbval % 生成的榜样比其邻居差  %生成的榜样比 pbest 差     
             iwt=0.75 + 0.005 * tan(pi * (rand(1, D) - 0.5));
             vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(pbest(j,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  （1）
           
        end
            
        vel(j,:)=(vel(j,:)>mv(1,:)).*mv(1,:)+(vel(j,:)<=mv(1,:)).*vel(j,:); 
        vel(j,:)=(vel(j,:)<(-mv(1,:))).*(-mv(1,:))+(vel(j,:)>=(-mv(1,:))).*vel(j,:);
        
        % 更新个体位置
        pos(j,:)= pos(j,:)+vel(j,:); 
        pos(j,:)=(pos(j,:)>VRmax(1,:)).*VRmax(1,:)+(pos(j,:)<=VRmax(1,:)).*pos(j,:); 
        pos(j,:)=(pos(j,:)<VRmin(1,:)).*(VRmin(1,:))+(pos(j,:)>=VRmin(1,:)).*pos(j,:);
        
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
        fitcount=fitcount+1;
        
        if e(j,1)<pbestval(j,1)
           pbestval(j,1)=e(j,1);
           pbest(j,:)=pos(j,:);
           if pbestval(j,1)<gbestval
               gbestval=pbestval(j,1);
               gbest=pbest(j,:);
           end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----        
    old2 = new2; 
    new2 = fitcount;
    de_dd_size=size(Q_df_dd_pos,1);
    center = mean(Q_df_dd_pos);
    div_i = 0;
    for di=1:de_dd_size    
        div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
    end
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
        ceased=0;
    else
        ceased=ceased+1;
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
%     end
 
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
    
% % diversity(Max_FES)=diversity(new2);
% xlabel('fes');
% ylabel('diversity');
% set(gca, 'Fontname', 'Times New Roman','FontSize',9);
% % x=1:round(Max_FES/100):Max_FES;
% x=1:round(6000/100):6000;
% % x(101)=Max_FES;
% % diversity(Max_FES) = gbestval;
% hold on;
% plot(x,diversity(x),'-r','MarkerFaceColor','r','MarkerSize',5);  
 
end

%%%%%%%%%%%%%%%%%% 更新档案  %%%%%%%%%%%%%%%%%%%%%%
function [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, max_size, obj_x, obj_val)
    archiev_size=length(archieve_val);
    if archiev_size<max_size  % 档案容量未达到上限值，直接放入 %
        archiev_size=archiev_size+1;
        archieve_val(archiev_size,1)=obj_val;
        archieve_pos(archiev_size,:)=obj_x;   
    else  % 档案容量已达到上限值，替换最差的个体 %
        [tmpB, tmpIX]=max(archieve_val);
        archiev_size=length(archieve_val);
%         if tmpB>obj_val  % obj_x比当前档案中最差个体要优时，替换该个体
        archieve_val(tmpIX,1)=obj_val;
        archieve_pos(tmpIX,:)=obj_x;
%         end

        % 档案容量已达到上限值，在档案中任选2个，替换其中较差的个体 %
%         rp=randperm(archiev_size);
%         if archieve_val(rp(1))<archieve_val(rp(2))% && archieve_val(rp(2))>obj_val % 更改了，注意！！！！！！        
%             archieve_val(rp(2),1)=obj_val;
%             archieve_pos(rp(2),:)=obj_x;
%                         
%         else%if archieve_val(rp(2))<archieve_val(rp(1)) && archieve_val(rp(1))>obj_val %  % 更改了，注意！！！！！！
%                 archieve_val(rp(1),1)=obj_val;
%                 archieve_pos(rp(1),:)=obj_x;            
%         end

    end
end

%%%%%%%%%%%%%%%%%% 重新利用档案更新个体的 pbest  %%%%%%%%%%%%%%%%%%%%%%
function [pbest, pbestval, vel]=Reuse_Exemplars(pbest, pbestval, vel, archieve_pos, archieve_val)
    archiev_size=length(archieve_val); 
    [ps,D]=size(pbest);
   
    [Val_p, Index_p]=sort(pbestval,'descend');
    i_p=1;
    
    index=1;
    rp=randperm(archiev_size);
    for i=1:floor(archiev_size/2)
        if archieve_val(rp(2*i-1))<archieve_val(rp(2*i))
           tmpval(index)=archieve_val(rp(2*i-1)); tmppos(index,:)=archieve_pos(rp(2*i-1),:);
           index=index+1;
        else
           tmpval(index)=archieve_val(rp(2*i)); tmppos(index,:)=archieve_pos(rp(2*i),:); 
           index=index+1;
        end
    end
    
    maxsize=length(tmpval);
    index=1;
    while i_p<=ps && index<=maxsize
        if Val_p(i_p,1)>tmpval(index) && rand<0.5
           pbestval(Index_p(i_p))=tmpval(index);
           pbest(Index_p(i_p),:)=tmppos(index,:);  

           vel(Index_p(i_p),:)=zeros(1,D);    
           i_p=i_p+1;
        end        
        index=index+1;
    end

end


%%%%%%%%%%%%%%%%%% 为子种群每个个体选择合适的学习对象  %%%%%%%%%%%%%%%%%%%%%%
function [obj_x,prob_f,prob_df,prob_df_dd]=select(i,ps,elites,Queue_f,pos_f,prob_f, Queue_df,pos_df,prob_df, Queue_df_dd, pos_df_dd, prob_df_dd,  Xmax,Xmin)
    pelite=0.20; % 精英个体的比例
    pcommon=0.40;% 平凡个体的比例
    pworst=0.40; % 劣质个体的比例

    [obj_x, prob_f, prob_df]=Generate(Queue_f,prob_f,pos_f,  Queue_df_dd,prob_df_dd,pos_df_dd, Xmax,Xmin);  % 
%     if i<=ceil(ps*pelite)  % 精英个体
%         [obj_x, prob_f, prob_df]=Generate(Queue_f,prob_f,pos_f,  Queue_df,prob_df,pos_df, Xmax,Xmin);  % A: 优质个体选择学习目标个体
% 
%     elseif i<=ceil(ps*(pelite+pcommon))  % 平凡个体       
%         [obj_x, prob_f, prob_df_dd]=Generate(Queue_f,prob_f,pos_f,  Queue_df_dd,prob_df_dd,pos_df_dd, Xmax,Xmin);  % B: 平凡个体选择学习目标个体
%  
%     else  % 劣质个体
%         [obj_x, prob_df_dd, prob_df]=Generate(Queue_df_dd,prob_df_dd,pos_df_dd, Queue_df,prob_df,pos_df,  Xmax,Xmin);  % C: 劣质个体选择学习目标个体
%     end
 
end

%%%%%%%%%%%%%%%%%%% 为某个体选择学习对象 %%%%%%%%%%%%%%%%%%%%%%%%
function [exemplar, prob_1, prob_2]=Generate(Queue_1,prob_1,pos_1, Queue_2,prob_2,pos_2,   Xmax,Xmin)
    N=size(Queue_1,1);
    D=size(pos_1,2);%返回pos_1的列数   
    pc=0.50; % 杂交概率
    pm=0.02; % 变异概率  
%    pm=0.02 + 0.005 * tan(pi * (rand(1, D) - 0.5));
    
    % 选择第1个个体时的概率 --------------------------
    p_1=prob_1;
    for i=2:N
       p_1(i)=p_1(i)+p_1(i-1);  % 计算第1个队列中各个体的候选概率
    end
    
    % 选择第2个个体时的概率 ---------------------------
    p_2=prob_2;
    for i=2:N
       p_2(i)=p_2(i)+p_2(i-1);  % 计算第2个队列中各个体的候选概率
    end
    
    
    %%
    for j=1:D    % 逐维选择学习榜样 
        for i=1:N
           if rand<=p_1(i)
               index=i;
               x_1=pos_1(index,:);
               break;
           end
        end    

        for i=1:N
           if rand<=p_2(i)
               index=i;
               x_2=pos_2(index,:);
               break;
           end
        end

        % 生成新个体 -------------------------------
        rd=rand;
        rd=rd<pc;   % 有此语句为几何杂交，无此语句则为算术杂交
        exemplar(1,j)=rd.*x_1(1,j)+(1-rd).*x_2(1,j);   % 利用2个个体算术杂交生成新的个体
    end
    
     %%
%       % 整体选择学习榜样 
%     for i=1:N
%        if rand<=p_1(i)
%            index=i;
%            x_1=pos_1(index,:);
%            break;
%        end
%     end    
% 
%     for i=1:N
%        if rand<=p_2(i)
%            index=i;
%            x_2=pos_2(index,:);
%            break;
%        end
%     end
% 
%     % 生成新个体 -------------------------------
%     rd=rand(1,D);
%     rd=rd<0.5;   % 有此语句为几何杂交，无此语句则为算术杂交
%     exemplar(1,:)=(1-rd).*x_1(1,:)+rd.*x_2(1,:);   % 利用2个个体算术杂交生成新的个体

    %% 变异
    for j=1:D
       if rand<pm
           exemplar(1,j)=Xmin(1,j)+(Xmax(1,j)-Xmin(1,j)).*rand;  
       end
    end
    
%     exemplar=pos_1(1,:);  %%%%%%%%%%%%%%%%%%%%%%%%%%% 词句为测试用！！！！
end

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
