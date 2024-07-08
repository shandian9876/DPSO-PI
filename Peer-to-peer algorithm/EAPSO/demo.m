clc
clear all
% mex cec13_func.cpp -DWINDOWS
nPop=50;
nVar=30;
MaxIt=150000;
% fhd=@sphere;
fhd=str2func('cec17_func');
VarMin=-100.*ones(1,nVar);
VarMax=100.*ones(1,nVar);
for i=1:nPop
    X(i,:)=VarMin+(VarMax-VarMin).*rand(1,nVar);
end
% [BestCost,BestValue] =EAPSO(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X);


runs=1; % 设置实验的运行次数为1
for i=        1:1
    func_num=i;
    for j=1:runs

        [BestCost,BestValue,diversity,everyfit] =EAPSO(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X,func_num);
        xbest(i,:)=BestCost;
        fbest(j,i)=BestValue;
        
        fprintf('第 %d 次运行的最优结果为：%1.4e\n',j,BestValue);
    end

    f_mean(func_num*2,:)=mean(fbest(:,func_num));
    f_mean(func_num*2+1,:)=std(fbest(:,func_num));
    
    fprintf('\nFunction F%d :\nAvg. fitness = %1.2e(%1.2e)\n\n',func_num,f_mean(func_num*2,:),f_mean(func_num*2+1,:));    
    fprintf(' -------------------------------------------------- \n');

end

xlabel('iteration');
ylabel('diversity');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'-x','color','m','MarkerFaceColor','m','MarkerSize',3,'LineWidth', 0.5);  
legend('DPSO-PI','EOPSO','TAPSO','EAPSO'); 


% plot(1:MaxIt/nPop*2,BestCost,'r')
% xlabel('The number of function evaluations','Fontname','Times New Roma','fontsize',15','FontWeight','bold');
% ylabel('Fitness value','Fontname','Times New Roma','fontsize',15,'FontWeight','bold');
