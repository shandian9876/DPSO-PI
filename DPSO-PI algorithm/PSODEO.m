function [Gbest_val,everyfit,diversity,CC]=PSODEO(fhd,MaxFes,Maxiter,N,D,Xmax,Xmin,optima_value,varargin)
% Knowledge-Based Systems期刊
count=zeros(N,1);
rand('state',sum(100*clock));
Vmax = Xmax/2;
Vmin = -Vmax;
fitness=zeros(N,1);
X = rand(N,D).*(Xmax-Xmin)+Xmin;
V = rand(N,D).*(Vmax-Vmin)+Vmin;
for i =1:N
    fitness(i,:) = feval(fhd,X(i,:)',varargin{:}); %feval(fhd,X(i,:)',varargin{:});
end
FES = N;
% Update pbest and gbest
[best_fitness,indFitness] = sort(fitness,'ascend'); 
Gbest = X(indFitness(1),:);               
Pbest = X;                            
Pbest_val = fitness;                
Gbest_val = best_fitness(1); 

everyfit=zeros(1,Maxiter);
diversity=zeros(1,Maxiter);
exemplar=zeros(N,D);
CC=zeros(1,Maxiter);
  %% 计算每个个体的 Novelty
% Novelty_type = 1;  
% num_neighbors=2;
% Novelty = Computing_Novelty(N, D, Pbest, Gbest, num_neighbors, Novelty_type); 
% [~, indNovelty] = sort(Novelty, 'descend'); 

%%
for iter = 1:Maxiter
 
    diversity(iter)=sum(sqrt(sum((X-mean(X)).^2,2)))/N;
    %=================== Find many solutions with higher fitness ==========
    p_ratio_fit = 0.6*(1./(1 + exp(0.0025*(iter-Maxiter/2))))+0.02; %[0.6 0.05]
%     p_ratio_fit=0.9*(Maxiter-iter)/Maxiter+0.02;
    pNP = p_ratio_fit.*rand([1 N]); 
    randindex1 = ceil(N * pNP); % select from [1, 2, 3, ..., pNP]
    randindex1 = max(1, randindex1); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    IX1 = indFitness(randindex1);
     pbest_fit = Pbest(IX1, :);
%     pbest_fit = (Pbest(IX1, :)+Gbest)/2;
  
        exemplar=DE_based_exemplar(N,D,Pbest,Xmax,Xmin);
        
    
         %=================== Find many solutions with higher novelty ==========
%     pNP2 = p_ratio_fit.*rand([1 N]);
%     randindex2 = ceil(N * pNP2); % select from [1, 2, 3, ..., pNP]
%     randindex2 = max(1, randindex2); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
%     IX2 = indNovelty(randindex2);
%     pbest_novelty = Pbest(IX2, :); % randomly choose one of the top 100p% solutions
    
%         [~,b]=sort(Pbest_val);
%         WW=(Pbest_val(b(10))-Pbest_val(b(1:10))+0.0001)/(Pbest_val(b(10))-Pbest_val(b(1))+0.0015);
%         PPbest=Pbest(b(1:10),:);
%         Mbest=sum(PPbest.*(WW/(sum(sum(WW)))));
        w=0.8-0.5*iter/Maxiter;%[0.9 0.4]
        c1=2.5-2*iter/Maxiter;%[2.5 0.5]
        c2=0.5+2*iter/Maxiter;%[0.5 2.5]
        cc(1) = 2*(1/(1 + exp(0.002*(iter-Maxiter/2)))); %单减[1.4 0.2]
        cc(2) = 0.5+ cc(1); %cc(2);% 单减[]
        cc(3) = 2.5-cc(1);%[]
 

        V=w.*V+c1*rand([N 1]).*(exemplar-X)+c2*rand([N 1]).*(pbest_fit-X);
        %+cc(2)*rand(N,D).*(Gbest-X);
%         V=w.*V+c1*rand([N 1]).*(pbest_novelty-X)+c2*rand([N 1]).*(pbest_fit-X);
       
        V=(V>Vmax).*Vmax+(V<=Vmax).*V;
        V=(V<Vmin).*Vmin+(V>=Vmin).*V;   
        X = X + V;    
        X=(X>Xmax).*Xmax+(X<=Xmax).*X;
        X=(X<Xmin).*Xmin+(X>=Xmin).*X;
        
        fitness = feval(fhd,X',varargin{:})';
        FES = FES + N;
      
        %=================== PSO: Update pbest and gbest ======================
    improved=(Pbest_val>fitness);     % individuals are improved or not
    temp=repmat(improved,1,D);
    Pbest=temp.*X+(1-temp).*Pbest;
    Pbest_val=improved.*fitness+(1-improved).*Pbest_val;      % update the pbest
    
    [val,indFitness]=sort(Pbest_val,'ascend');
    Gbest_val=val(1);
    Gbest=Pbest(indFitness(1),:);              % update the gbest

    
    everyfit(iter)=Gbest_val;
    if Gbest_val==optima_value
        Gbest_val=optima_value;
        for kkk=iter:Maxiter
            everyfit(kkk)=Gbest_val;
        end
        break;
    end
    
      %% 计算每个个体的 Novelty
%     Novelty = Computing_Novelty(N, D, Pbest, Gbest, num_neighbors, Novelty_type); 
%     [~, indNovelty] = sort(Novelty, 'descend'); 

     %%
    % % 跳出策略
    if  iter < 0.8*Maxiter 
            V= (abs(V)>10^(-4)).* V+(1-(abs(V)>10^(-4))).*normrnd(0,1,[1,D]);
    end
    
    flag = 1-improved;
    count = flag + flag.*count;
    fiter = 10^floor(log10(Maxiter));
    K = 60 + 40*tanh((iter/fiter)-(Maxiter/fiter)/2); %[20，100] 

    flog = (count >= K);
%     if(any(flog>=1))
    if(1==2)
        count = (1-flog).*count;
        p = rand(N,D).*(Xmax-Xmin)+Xmin;
        X = (1-flog).*X + flog.*p;
        
         fitness = feval(fhd,X',varargin{:})';
      
        %=================== PSO: Update pbest and gbest ======================
        improved=(Pbest_val>fitness);     % individuals are improved or not
        temp=repmat(improved,1,D);
        Pbest=temp.*X+(1-temp).*Pbest;
        Pbest_val=improved.*fitness+(1-improved).*Pbest_val;      % update the pbest
        
        [val,indFitness]=sort(Pbest_val,'ascend');
        Gbest_val=val(1);
        Gbest=Pbest(indFitness(1),:);              % update the gbest
    end

end
end
