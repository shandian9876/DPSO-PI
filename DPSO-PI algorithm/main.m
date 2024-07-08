clc
clear all
% mex cec17_func.cpp -DWINDOWS

D=30;
Xmin=-100;
Xmax=100;
pop_size=50;
MaxFES=10000*D;
Maxiter=ceil(MaxFES/pop_size);

% % 预定义的优化目标值
optima_value=[100, 200, 300, 400, 500,...
       600, 700, 800, 900, 1000,...
       1100,1200,1300,1400,1500,...
       1600,1700,1800,1900,2000,...
       2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000 ];
  
% 设置实验的运行次数为1 
runs=10; 
fhd=str2func('cec17_func');
for i=       1:30
    func_num=i;
    for j=1:runs  
        [Gbest_val,everyfit,diversity,CC]=PSODEO(fhd,MaxFES,Maxiter,pop_size,D,Xmax,Xmin,optima_value(func_num),func_num);

%         xbest(j,:)=gbest;
        fbest(j,func_num)=Gbest_val;
        fprintf('第 %d 次运行的最优结果为：%1.4e\n',j,Gbest_val);
    end
%     f_mean(i)=mean(fbest(i,:));
    f_mean(func_num,:)=mean(fbest(:,func_num));
%     f_mean(func_num*2,:)=std(fbest(:,func_num));
    fprintf('\nFunction F%d :\nAvg. fitness = %1.2e(%1.2e)\n\n',func_num,mean(fbest(:,func_num)),std(fbest(:,func_num)));    
    fprintf(' -------------------------------------------------- \n');
 end

%%convergence speed
xlabel('iteration');
ylabel('fitness');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'-*','color','r','MarkerFaceColor','r','MarkerSize',3,'LineWidth', 0.5);  
legend('DPSO-PI');

%%diversity
% xlabel('iteration');
% ylabel('diversity');
% set(gca, 'Fontname', 'Times New Roman','FontSize',9);
% x=1:round(6000/100):6000;
% hold on;
% plot(x,diversity(x),'-r','MarkerFaceColor','r','MarkerSize',8,'LineWidth', 1.2);  
% legend('DPSO-PI');




