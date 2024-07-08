function [Gbest_val,everyfit,diversity]= EOPSO(fhd,MaxFES,N,D,LB,UB,varargin)
% Gbest_val:全局最佳位置的适应度
% Pbest_val:个人最佳位置的适应度
% fitness_record:存储收敛曲线
% MaxFES Maximum:适应度评估次数
% swarm size:群体大小
% dimension:问题维度
% LB/UP:上/下搜索边界
% pos:粒子位置
% inw:惯性权重
% diversity是计算每次迭代中群体的多样性度量值，并将其存储到 diversity 数组中。用于论文中多样性分析；
% 通过跟踪多样性度量值的变化，可以了解群体在不同迭代中的多样性水平，以评估算法的收敛性
% 和搜索空间的探索能力。多样性度量值越大，表示群体的多样性越高，搜索空间的探索能力越强。
rand('state',sum(100*clock));


Maxiter=ceil(MaxFES/N);
everyfit=zeros(1,Maxiter);
velocityelocity_max=(UB-LB)/4;
velocityelocity_min=-velocityelocity_max;
pos=rand(N,D).*(UB-LB)+LB;
velocity=rand(N,D).*(velocityelocity_max-velocityelocity_min)+velocityelocity_min;
fitness=zeros(N,1);
for i=1:N
    fitness(i)=feval(fhd,pos(i,:)',varargin{:});
end
fitcount=N;
everyfit(1)=min(fitness);
[min_val,min_index]=min(fitness);
Gbest=pos(min_index,:);
Gbest_val=min_val;
Pbest=pos;
Pbest_val=fitness;
count=zeros(1,N);
c=1.49445;
iter=1;
diversity=zeros(1,Maxiter);
K1=0.5;
K2=0.4;
G=7; %如果粒子停止更新G次迭代，则进行跳出策略

while iter<=Maxiter %& fitcount<=MaxFES


    diversity(iter)=sum(sqrt(sum((pos-mean(pos)).^2,2)))/N;
    inw=0.9-0.5*(iter/Maxiter);
    iter=iter+1;

    elite_num=ceil(N*K1-iter/Maxiter*(K1-K2)*N);

    [~,index]=sort(fitness);
    OO=zeros(1,D);
    for i=1:elite_num
        if fitness(index(i))>0
            f(i)=1/(fitness(index(i))+1);
        else
            f(i)=1+abs(fitness(index(i)));
        end
    end

    elite=zeros(1,elite_num);
    elite(1:elite_num)=index(1:elite_num);
    ordinary=zeros(1,N-elite_num);
    ordinary(1:N-elite_num)=index(elite_num+1:N);

    for i=1:elite_num
        OO=OO+f(i)/sum(f)*Pbest(index(i),:);
    end

    for i=1:N

        ranknum=index(i);

        if ranknum<=elite_num
            velocity(i,:)=inw.*velocity(i,:)+c*rand(1,D).*(Pbest(i,:)-pos(i,:));
        else
            velocity(i,:)=inw.*velocity(i,:)+c*rand(1,D).*(OO-pos(i,:));
        end

        velocity(i,velocity(i,:)>velocityelocity_max) = velocityelocity_max;
        velocity(i,velocity(i,:)<velocityelocity_min) = velocityelocity_min;
        pos(i,:) = pos(i,:) + velocity(i,:);
        pos(i,pos(i,:)>UB) = UB;
        pos(i,pos(i,:)<LB) = LB;
        fitness(i)=feval(fhd,pos(i,:)',varargin{:});
        fitcount=fitcount+1;
%         everyfit(fitcount)=min(everyfit(fitcount-1),fitness(i));


        Pbest(i,:)=(fitness(i)<Pbest_val(i))*pos(i,:)+(fitness(i)>=Pbest_val(i))*Pbest(i,:);
        Pbest_val(i)=(fitness(i)<Pbest_val(i))*fitness(i)+(fitness(i)>=Pbest_val(i))*Pbest_val(i);
        count(i)=(fitness(i)>Pbest_val(i))+(fitness(i)>Pbest_val(i))*count(i);

        Gbest=(fitness(i)<Gbest_val)*pos(i,:)+(fitness(i)>=Gbest_val)*Gbest;
        Gbest_val=(fitness(i)<Gbest_val)*fitness(i)+(fitness(i)>=Gbest_val)*Gbest_val;

        everyfit(iter)=Gbest_val;%min(fitness);

        if count(i)>=G
            count(i)=0; 
            r1=(ranknum<=elite_num)*ordinary(randperm(N-elite_num,1))+(1-(ranknum<=elite_num))*elite(randperm(elite_num,1));
            RR=rand(1,D);
            Ppos=RR.*Pbest(i,:)+(1-RR).*Pbest(r1,:)+2*(rand(1,D)-0.5).*(Pbest(i,:)-Pbest(r1,:));
            Ppos_fit=feval(fhd,Ppos',varargin{:});
            fitcount=fitcount+1;
%             everyfit(fitcount)=min(everyfit(fitcount-1),Ppos_fit);

            pos(i,:)=Ppos*(Ppos_fit<fitness(i))+pos(i,:)*(1-(Ppos_fit<fitness(i)));
            fitness(i)=Ppos_fit*(Ppos_fit<fitness(i))+fitness(i,:)*(1-(Ppos_fit<fitness(i)));
            Pbest(i,:)=Ppos*(Ppos_fit<Pbest_val(i))+Pbest(i,:)*(1-(Ppos_fit<Pbest_val(i)));
            Pbest_val(i)=Ppos_fit*(Ppos_fit<Pbest_val(i))+Pbest_val(i,:)*(1-(Ppos_fit<Pbest_val(i)));
        end
        if  iter<=0.9*Maxiter
            velocity(i,:)= (abs(velocity(i,:))>10^(-4)).* velocity(i,:)+(1-(abs(velocity(i,:))>10^(-4))).*normrnd(0,1,[1,D]);
        end
    end

%     if fitcount>=MaxFES
%         break;
%     end
%     if (iter==Maxiter) & (fitcount<MaxFES)
%         iter=iter-1;
%     end
end
% everyfit=everyfit(1:MaxFES);
% Gbest_val=everyfit(1,MaxFES);

end

