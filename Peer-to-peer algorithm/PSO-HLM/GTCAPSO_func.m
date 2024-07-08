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
%%%%%%%%%%%%%%%%%% ��������� --- ������ ----
% old2 = 1; new2 = fitcount;  %%%------��¼���۴����ĸı�
% de_dd_size=size(Q_df_dd_pos,1);
% center = mean(Q_df_dd_pos);     
% div_i = 0;      %%%
% for di=1:de_dd_size    
%     div_i = div_i + sum((Q_df_dd_pos(di,:)-center).^2);  
% end             %%%������ά�ȵ�λ�õ�ֵ���
% diversity(old2:new2) = div_i/(de_dd_size*(VRmax(1,1)-VRmin(1,1))) ; %%% �� ������   d��j�� =�� X ij - ƽ����X j����1��m���ۼӺ�
% %%%%%%%%%%%%%%%%%%

old_e=e; old_pos=pos;        %%% �����ڵ�λ�ñ�Ϊ��һ����λ�ã��Ա�������� λ�ñ仯
old_pbestval=pbestval; old_pbest = pbest;       %%% �����ڵ���Ӧ�ȱ�Ϊ��һ������Ӧ�ȣ��Ա�������� ��Ӧ�ȱ仯

recorded = 0;  % �ﵽ����ʱ��¼�����Ϣ
suc = 0;        %%%-----
suc_fes = 0;    %%%-----

%%%% ����Ϊ��ͼ��
old = 1;
new = fitcount;
yyy(old:new) = gbestval;
old = new;
iter=1;
ceased=0;

archieve_val=[];%%%��ʼ���������� ��Ӧ�� ����
archieve_pos=[];%%%��ʼ���������� λ�� ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jg(fitcount)=gbestval;%���������
% while fitcount<Max_FES
while iter<me
    diversity(iter)=sum(sqrt(sum((pos-mean(pos)).^2,2)))/ps;
    everyfit(iter)=gbestval;
    iter=iter+1;  %%%  iter ��+1����Ϊ�ڸտ�ʼ���Ѿ����¹�һ��  ���Դӵڶ��ο�ʼ 
%%%-------------- ������������-------------
    for k=1:ps
        % ѡ��ѧϰĿ�����          
        [obj_x(k,:),prob_f,prob_df,prob_df_dd]=select(iter,me,k,ps,sizeof_Af, Q_f,Q_f_pos,prob_f, Q_df,Q_df_pos,prob_df, Q_df_dd,Q_df_dd_pos,prob_df_dd, VRmax,VRmin);       
%        ��ʱ obj_x �ɽ������ۣ������ݽ���������뵽��Ⱥ��ȥ
        obj_val(k)=feval(fhd,obj_x(k,:)',varargin{:});  %%%������в�������
        if iter>1/3*me
            obj_val1(k)=obj_val(k)*abs(randn);
        else
            obj_val1(k)=obj_val(k);
        end
        fitcount=fitcount+1;    %%%������֮�� ����+1������¼++++++++++++++++++++++++
         jg(fitcount)=gbestval;%�������
    end 
%%   
    [~, IX]=sort(obj_val1); 
       for k=1:ps
           j=IX(k); %%% ��IX��1����ʼѡ��   ���� gbest �𽥽��� ��ʼѡ��
        if obj_val(j)<=gbestval  % ���ɵİ����� gbest ��
            [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x(j,:), obj_val(j));
            iwt=0.45 + 0.005 * tan(pi * (rand(1, D) - 0.5));        %%%���㱣��ϵ��
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(3).*rand(1,D).*(obj_x(j,:)-pos(j,:));   % ��3��         %%% ���� ��ɢ�ͽ����ٶȸ���
            gbest(1,:)=obj_x(j,:);
            gbestval=obj_val(j);
            
         elseif obj_val1(j)<=gbestval %nbval% ���ɵİ��������ھӺ�  %���ɵİ����� pbest �� 
               iwt=0.45 + 0.005 * tan(pi * (rand(1, D) - 0.5));        %%%���㱣��ϵ��
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(3).*rand(1,D).*(obj_x(j,:)-pos(j,:))+0.001*abs(randn);   % ��3��         %%% ���� ��ɢ�ͽ����ٶȸ���
            
        elseif obj_val1(j)<=pbestval(j)  %nbval% ���ɵİ��������ھӺ�  %���ɵİ����� pbest �� 
           [archieve_val, archieve_pos]=add2archieve(archieve_val, archieve_pos, archive_maxsize, obj_x(j,:), obj_val(j));
            iwt=0.7298 + 0.005 * tan(pi * (rand(1, D) - 0.5));       %%%���㱣��ϵ��
            vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(obj_x(j,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  ��2��%%% ���� �����ͽ����ٶȸ���
           
        elseif obj_val1(j)>pbestval(j)  % nbval % ���ɵİ��������ھӲ�  %���ɵİ����� pbest ��     
             iwt=0.7298+ 0.005 * tan(pi * (rand(1, D) - 0.5));         %%%���㱣��ϵ��
             vel(j,:)=iwt(1,:).*vel(j,:)+cc(1).*rand(1,D).*(pbest(j,:)-pos(j,:))+cc(2).*rand(1,D).*(gbest(1,:)-pos(j,:));  %  ��1��%%% ���� �º��ͽ����ٶȸ���    
        end
            
        vel(j,:)=(vel(j,:)>mv(1,:)).*mv(1,:)+(vel(j,:)<=mv(1,:)).*vel(j,:);     %%% ���Ƹ����ٶ�
        vel(j,:)=(vel(j,:)<(-mv(1,:))).*(-mv(1,:))+(vel(j,:)>=(-mv(1,:))).*vel(j,:);  %%%���Ƹ����ٶ�
        
        % ���¸���λ��
        pos(j,:)= pos(j,:)+vel(j,:);  %%%���¸���λ��
        pos(j,:)=(pos(j,:)>VRmax(1,:)).*(VRmax(1,:))+(pos(j,:)<=VRmax(1,:)).*pos(j,:);   %%% ���Ƹ���λ��
        pos(j,:)=(pos(j,:)<VRmin(1,:)).*(VRmin(1,:))+(pos(j,:)>=VRmin(1,:)).*pos(j,:);   %%%���Ƹ���λ��
        
         % Խ�紦��
        if rand<0.5
            pos(j,:)=(pos(j,:)>VRmax(1,:)).*VRmax(1,:)+(pos(j,:)<=VRmax(1,:)).*pos(j,:); 
            pos(j,:)=(pos(j,:)<VRmin(1,:)).*(VRmin(1,:))+(pos(j,:)>=VRmin(1,:)).*pos(j,:);
        else
            pos(j,:)=((pos(j,:)>=VRmin(1,:))&(pos(j,:)<=VRmax(1,:))).*pos(j,:)...
                +(pos(j,:)<VRmin(1,:)).*(VRmin(1,:)+0.1.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D))...
                +(pos(j,:)>VRmax(1,:)).*(VRmax(1,:)-0.1.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D));           
        end
        %���۸���
        e(j,1)=feval(fhd,pos(j,:)',varargin{:}); 
%         e(j,1)=feval(fhd,pos(j,:)',varargin{:})-fbias(func_num);
        fitcount=fitcount+1;        %%%+++++++++++++++++++++++++++++++++
         jg(fitcount)=gbestval;%�������
         
        if e(j,1)<pbestval(j,1)
           pbestval(j,1)=e(j,1);
           pbest(j,:)=pos(j,:);
           if pbestval(j,1)<gbestval
               gbestval=pbestval(j,1);
               gbest=pbest(j,:);
           end
        end
        
       end
       
      jg(fitcount)=gbestval;%���������
      
    %%%%%%%%%%%%%%%%%% ��������� --- ������ ----        
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
%    % �� archieve �еĸ����� pbest �ȽϺ�̰������
  [pbest, pbestval, vel]=Reuse_Exemplars(pbest, pbestval, vel, archieve_pos, archieve_val);
   
   % ���������ǽ���Ⱥ������ʷ�����еĸ�����ʷ���Ž���뵵����
   % �п��ܻὫ��Ӣ����Ķ�θ�����ʷ���Ž����õ������Ƿ������ν���
   % ���⣬����� oldpbest �ƺ�Ҳû����ÿ���н��и��£���Ҫ�޸Ĳ��ִ���
   % ��ע�͵�����һ�䣻Ȼ���� Update_Q ����������� old_pbestval = pbestval �ȼ������
   
   
%%
    % ������Ⱥ��ʷ���Ž�
    [~, Index]=min(pbestval);
    if pbestval(Index)<gbestval 
        gbest=pbest(Index,:);
        gbestval=pbestval(Index);
%         gbestrep=repmat(gbest,ps,1);%update the gbest   
        ceased=0;   %%%-----------
    else
        ceased=ceased+1;        %%%----------------
    end  
   
      
    % ������Ӧֵ�仯��Ϣ ------------   
     pbestval_df=old_pbestval - pbestval;  % �� pbest �ı仯Ϊָ�� 
%      pbestval_df=old_e - e;  % �� pos �ı仯Ϊָ��
    
%     [pbestval_dif IX]=sort(pbestval_dif,'descend');  % ��Ӧֵ��������
 
    % ���� ��Ӧֵ�仯/����仯 ��Ϣ
%     df=old_e-e; % ��Ӧֵ�仯
%     dd=sqrt(sum((pos(:,:)-old_pos(:,:)).*(pos(:,:)-old_pos(:,:)),2)); % λ�ñ仯
    
    df=old_pbestval - pbestval;   % ʹ�� pbestval 
    dd=sqrt(sum((pbest(:,:)-old_pbest(:,:)).*(pbest(:,:)-old_pbest(:,:)),2)); % pbestλ�ñ仯
    
    df_dd=df./dd; % ����Ч��
%     [tmp_df_dd IX]=sort(tmp_df_dd,'descend');  % ����Ч�ʽ�������
%     tmppos_df_dd=pos(IX(:),:);          % ��Ӧ����Ҳ�ź���
  
    
    % ���� Q_f,Q_df,Q_df_dd ��ÿ�������и��������Ȩ
    [Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd]=...
        Update_Q(Q_f,Q_f_pos,prob_f,  Q_df,Q_df_pos,prob_df,  Q_df_dd,Q_df_dd_pos,prob_df_dd, pbestval,pbest,  pbestval_df, df_dd,pos);
      
    % ����ӵ�
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
%%%%%%%%%%%%%%%%%% ��������� --- ������ ----        
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

