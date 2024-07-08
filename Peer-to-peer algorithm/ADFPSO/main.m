clear all
clc

global fbias
warning('off')

D=30;
Xmin=-100;
Xmax=100;

pop_size=50;%2^(4+floor(log2(sqrt(D))));
fes_max=10000*D;
iter_max=ceil(fes_max/pop_size);

runtimes=1;
fhd=str2func('cec17_func');


jingdu=0;

funset=[1:30];
'ADFPSO-CEC2017'  
% for fun=28:length(funset)
for fun=        1:1
    func_num=funset(fun);
    suc_times = 0; 
    
    fesusage=0;
    count=0;
    
    for runs=1:runtimes
        suc = 0;
        suc_fes = 0;   % added by us
        tic;
    
        [gbest,gbestval,FES,diversity,everyfit]= ADFPSO_func(jingdu,func_num,fhd,D,pop_size,iter_max,fes_max,Xmin,Xmax,func_num);

        t=toc;
        time_usage(runs,fun)=t;
        
        xbest(runs,:)=gbest;
        fbest(runs,func_num)=gbestval;
        fprintf('第 %d 次运行的最优结果为：%1.4e\n',runs,gbestval);
        suc_times = suc_times + suc;
        if suc == 1  % 当达到设定精度时才统计其耗时           
            fesusage = fesusage + suc_fes;   % 达到精度时的 fes            
        end        
    end

    f_mean(func_num*2,:)=mean(fbest(:,func_num));
    f_mean(func_num*2+1,:)=std(fbest(:,func_num));
    fprintf('\nFunction F%d :\nAvg. fitness = %1.2e(%1.2e)\n\n',func_num,mean(fbest(:,func_num)),std(fbest(:,func_num)));    
    fprintf(' -------------------------------------------------- \n');
    
end

xlabel('iteration');
ylabel('log10(fitness)');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'->','color','c','MarkerFaceColor','c','MarkerSize',3,'LineWidth', 0.5);  
legend('DPSO-PI','EOPSO','TAPSO','EAPSO','PSOHLM','ADFPSO'); 
%   ylim([0 40000])

