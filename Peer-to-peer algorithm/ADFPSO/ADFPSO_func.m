function [gbest,gbestval,FES,diversity,everyfit]= ADFPSO_func(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,MaxFes,VRmin,VRmax,varargin)
                                              
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)

% fbias=[100, 200, 300, 400, 500,...
%        600, 700, 800, 900, 1000,...
%        1100,1200,1300,1400,1500,...
%        1600,1700,1800,1900,2000,...
%        2100,2200,2300,2400,2500,...
%        2600,2700,2800,2900,3000 ];
fbias=[0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0 ];
n=Dimension;
D=n;

Xmin=VRmin; 
Xmax=VRmax;
% lu 的最后3列分别为 iwt, c1, c2 的上下界。因为使用 self-adaptive 策略
% lu=[Xmin* ones(1, D), iwt_min, c1_min, c2_min; Xmax* ones(1, D), iwt_max, c1_max, c2_max];
lu=[Xmin* ones(1, D); Xmax* ones(1, D)];

NP=Particle_Number;
Ninit=NP;

% fhd=str2func('cec13_func');

totalTime = 1;
bestresults=[];

% err = accuracy(problem);
err = jingdu;

% Record the best results
outcome = []; 
rand('seed', sum(100 * clock));

% Initialize the main population
FES = 0;
% pop = repmat(lu(1, :), NP, 1) + rand(NP, D+3) .* (repmat(lu(2, :) - lu(1, :), NP, 1));
pop = repmat(lu(1, :), NP, 1) + rand(NP, D) .* (repmat(lu(2, :) - lu(1, :), NP, 1));
popold = pop;

% Initialize the velocity
mv = 0.1*(lu(2,:) - lu(1,:));
Vmin=repmat(-mv,NP,1);
Vmax=-Vmin;
% vel=Vmin+2.*Vmax.*rand(NP,D+3);
vel=Vmin+2.*Vmax.*rand(NP,D);

% Set parameters
iwt=0.9-(1:Max_Gen)*(0.5/Max_Gen);
cc=[2.0, 2.0]; 

% valParents = benchmark_func(popold, problem, o, A, M, a, alpha, b);        
e=(feval(fhd,pop',varargin{:})-fbias(func_num))';
FES = FES + NP;

% Update pbest and gbest
pbest=pop;
pbestval=e; %initialize the pbest and the pbest's fitness value
[val,indFitness]=sort(pbestval,'ascend');
gbestval=val(1);
gbest=pbest(indFitness(1),:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,NP,1);

%%%% 以下为画图用    
old = 1;
new = 1;
yyy(old:new) = gbestval;     
%%%%

% %%%%%%%%%%%%%%%%%% 计算多样性 --- 测试用 ----
old2 = 1; new2 = FES;
pop_size=size(pop,1);
center = mean(pop);
div_i = 0;
for di=1:pop_size    
    div_i = div_i + sum((pop(di,:)-center).^2);  
end
% diversity(old2:new2) = div_i/(pop_size*(VRmax(1,1)-VRmin(1,1))) ; 

%% 计算每个个体的 Novelty
num_neighbors = 2;%ceil(Particle_Number*0.02);%   % 每个个体有 3 个邻居

% 1: 与最近的若干邻居间的距离表征 novelty; 2: 与 gbest 之间的距离; 
% 3: 与 pbest 中心点的距离; 4: 与 gbest 及 pbest 的中心点的距离之和
Novelty_type = 1;  
Novelty = Computing_Novelty(NP, D, pbest, gbest, num_neighbors, lu, Novelty_type); 
[~, indNovelty] = sort(Novelty, 'descend');

%%
Afactor = 2.6;
archive.NP = ceil(Afactor * NP); % the maximum size of the archive
archive.pop = zeros(0, n); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions

distance = zeros(0,1);
Ndelete = ceil(Ninit/(50));%ceil(Ninit/50);
Nmin = ceil(D);
%Nmin = ceil(Ninit/20);
oldlu = lu;
%%
recorded = 0;  % 达到精度时记录相关信息
suc = 0;
suc_fes = 0;
k1=1; k2=1;

g=0;
weight=1.0; % 当 pbest 进化时，其Novelty 乘以权重系数 weight
while g <  6001 %&& g < iter_max %n * 10000 %& min(fit)>error_value(problem)
    g = g + 1;   
    everyfit(g)=gbestval;
    diversity(g)=sum(sqrt(sum((pop-mean(pop)).^2,2)))/NP;
    
    p_ratio_fit = 0.7*(1/(1+exp(0.001*(FES-(MaxFes/2))/D)))+0.1;  % sigmoid 函数: [0.8, 0.1]
    p_ratio_novelty = p_ratio_fit;  % 0.9-p_ratio_fit;  %   
        
    %=================== Find many solutions with higher fitness ==========
    pNP = p_ratio_fit.*rand([1 NP]); 
    randindex1 = ceil(NP * pNP); % select from [1, 2, 3, ..., pNP]
    randindex1 = max(1, randindex1); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    IX1 = indFitness(randindex1);
    pbest_fit = pbest(IX1, :); % randomly choose one of the top 100p% solutions
   
    %=================== Find many solutions with higher novelty ==========
    pNP2 = p_ratio_novelty.*rand([1 NP]);
    randindex2 = ceil(NP * pNP2); % select from [1, 2, 3, ..., pNP]
    randindex2 = max(1, randindex2); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    IX2 = indNovelty(randindex2);
    pbest_novelty = pbest(IX2, :); % randomly choose one of the top 100p% solutions
    
    %=================== PSO: Update position and velocity ================
    iwt=0.9-(FES)*(0.5/MaxFes);    
    cc(2)=1.2*(1/(1+exp(0.0015*(FES-(MaxFes/2))/n)))+0.2;%1.5-(FES)*(1.2/MaxFes); %
    cc(1)=1.4-cc(2); %cc(2);% 
    aa=cc(1).*rand(NP,D).*(pbest_fit-pop)+cc(2).*rand(NP,D).*(pbest_novelty-pop);
    vel=iwt*vel+aa;
    
    %=================== PSO: Constrain position and velocity ================
    vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
    vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
    pop=pop+vel;  
    if rand>0.5
        pop=(pop>VRmax).*VRmax+(pop<=VRmax).*pop; 
        pop=(pop<VRmin).*VRmin+(pop>=VRmin).*pop;
    else
        pop=((pop>=VRmin)&(pop<=VRmax)).*pop ...
            +(pop<VRmin).*(VRmin+0.2.*(VRmax-VRmin).*rand(NP,D))+(pop>VRmax).*(VRmax-0.2.*(VRmax-VRmin).*rand(NP,D));
    end
    
    e=(feval(fhd,pop',varargin{:})-fbias(func_num))';
    FES=FES+NP;
    
    %=================== PSO: Update pbest and gbest ======================
    improved=(pbestval>e);     % individuals are improved or not
    temp=repmat(improved,1,D);
    pbest=temp.*pop+(1-temp).*pbest;
    pbestval=improved.*e+(1-improved).*pbestval;      % update the pbest
    
    [val,indFitness]=sort(pbestval,'ascend');
    gbestval=val(1);
    gbest=pbest(indFitness(1),:);              % update the gbest
    gbestrep=repmat(gbest,NP,1);

   %% 计算每个个体的 Novelty

    % 1: 与最近的若干邻居间的距离表征 novelty; 2: 与 gbest 之间的距离; 
    % 3: 与 pbest 中心点的距离; 4: 与 gbest 及 pbest 的中心点的距离之和
    Novelty_type = 1;  
    Novelty = Computing_Novelty(NP, D, pbest, gbest, num_neighbors, lu, Novelty_type); 
  
    [~, indNovelty] = sort(Novelty, 'descend');
  
    %%
    %%%%%%%%%%%%%%%    
    if gbestval <= jingdu && recorded == 0
        recorded = 1;
        suc = 1;
        suc_fes = FES;
%         break;
    end

end
  
end


%==========================================================================
function Novelty = Computing_Novelty(NP, D, pbest, gbest, num_neighbors, lu, type)
    switch type
    case 1
    % 个体的 pbest 与最临近的 num_neighbors 个个体 pbest 的距离来表征 Novelty
         for i = 1:NP
            for j = 1:NP
               distance(1,j) = norm(pbest(i,:)-pbest(j,:));
            end        

            distance(distance==0)=[]; % 删掉 0 元素，即不考虑个体与其自身的距离
            if ~isempty(distance)
                [Novel_val, Novel_IX] = sort(distance,'ascend');

                if length(Novel_IX) == 1   % 当种群已经收敛时:对 pbest 随机扰动
%                     Novelty(i,1) = mean(Novel_val(Novel_IX(1)));
%                     sigma = (lu(2)-lu(1))/200;
%                     pbest(i,:) = pbest(i,:) + normrnd(0,sigma,1,D);

                else
                    Novelty(i,1) = mean(Novel_val(Novel_IX(1:min(num_neighbors,length(Novel_IX)))));
                end               
            else
                Novelty(i,1) = 0;
            end
         end
     
         if max(Novelty) > 0
            Novelty = Novelty/max(Novelty);     % 归一化，值越大表明越新颖！
         else
            return;
         end
     
    case 2
         % 个体的 pbest 与gbest 的距离来表征 Novelty     
         for i = 1:NP        
            distance(i) = norm(pbest(i,:)-gbest(1,:));              
         end

         if sum(distance)==0  % 种群已收敛，对 pbest 进行扰动        
                sigma = (lu(2)-lu(1))/200;
                pbest(i,:) = pbest(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pbest(i,:)-gbest(1,:));              
                end
         end

         max_Novelty = max(distance);    
         Novelty = distance/max_Novelty;     % 归一化，值越大表明越新颖！   
            
     case 3
        % 个体的 pbest 与 pbest 的中心点的距离来表征 Novelty
         cbest = mean(pbest);
         for i = 1:NP        
            distance(i) = norm(pbest(i,:)-cbest(1,:));              
         end

         if sum(distance)==0  % 种群已收敛，对 pbest 进行扰动        
                sigma = (lu(2)-lu(1))/200;
                pbest(i,:) = pbest(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pbest(i,:)-gbest(1,:));              
                end
         end

         max_Novelty = max(distance);    
         Novelty = distance/max_Novelty;     % 归一化，值越大表明越新颖！  
         
      case 4
          % 个体的 pbest 与 gbest 及 pbest 的中心点的距离来表征 Novelty
          cbest = mean(pbest);
          for i = 1:NP        
            distance_1(i) = norm(pbest(i,:)-gbest(1,:));
            distance_2(i) = norm(pbest(i,:)-cbest(1,:));
            distance(i) = distance_1(i)+distance_1(i);
          end

          if sum(distance)==0  % 种群已收敛，对 pbest 进行扰动        
                sigma = (lu(2)-lu(1))/200;
                pbest(i,:) = pbest(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pbest(i,:)-gbest(1,:));              
                end
          end

          max_Novelty = max(distance);    
          Novelty = distance/max_Novelty;     % 归一化，值越大表明越新颖！
    end
end
