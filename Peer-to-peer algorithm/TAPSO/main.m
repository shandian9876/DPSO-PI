
clear all
clc

global fbias
warning('off')

D=30;
Xmin=-100;
Xmax=100;

pop_size=50;%2^(4+floor(log2(sqrt(D))));
N=pop_size;
fes_max=10000*D;
iter_max=ceil(fes_max/pop_size);



fhd=str2func('cec17_func');
% fhd=str2func('cec13_func');

fbias=[100, 200, 300, 400, 500,...
       600, 700, 800, 900, 1000,...
       1100,1200,1300,1400,1500,...
       1600,1700,1800,1900,2000,...
       2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000 ];

% fbias=[-1400, -1300, -1200, -1100, -1000,...
%        -900, -800, -700, -600, -500,...
%        -400,-300,-200,-100,100,...
%        200,300,400,500,600,...
%        700,800,900,1000,1100,...
%        1200,1300,1400 ];

jingdu=0;

funset=[1:30];
'ADFPSO-CEC2017'  

runtimes=1;
for fun=       1:1
%length(funset)
    func_num=funset(fun);
    suc_times = 0; 
    
    fesusage=0;
    count=0;
    jingdu=0;
    
    for runs=1:runtimes
        suc = 0;
        suc_fes = 0;   % added by us
        tic;
    
        [gbest,gbestval,fitcount,suc,FES,diversity,everyfit]= TAPSO_func(jingdu,fhd,iter_max,fes_max,N,D,-100,100,fun); 

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
ylabel('diversity');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'-o','color','b','MarkerFaceColor','b','MarkerSize',3,'LineWidth', 0.5);  

legend('DPSO-PI','EOPSO','TAPSO');




