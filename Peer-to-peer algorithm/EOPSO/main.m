clc
clear all
% mex cec17_func.cpp -DWINDOWS
func_num=1;
D=30;
Xmin=-100;
Xmax=100;
pop_size=50;
MaxFES=10000*D;
iter_max=ceil(MaxFES/pop_size);

runs=1; % 设置实验的运行次数为1
fhd=str2func('cec17_func');
for i=         1:1
    func_num=i;
    for j=1:runs
        [Gbest_val,everyfit,diversity] = EOPSO(fhd,MaxFES,pop_size,D,Xmin,Xmax,func_num);
%         xbest(j,:)=gbest;
        fbest(j,func_num)=Gbest_val;
        fprintf('第 %d 次运行的最优结果为：%1.4e\n',j,Gbest_val);
    end
%     f_mean(i)=mean(fbest(i,:));
    f_mean(func_num*2,:)=mean(fbest(:,func_num));
    f_mean(func_num*2+1,:)=std(fbest(:,func_num));
    fprintf('\nFunction F%d :\nAvg. fitness = %1.2e(%1.2e)\n\n',func_num,mean(fbest(:,func_num)),std(fbest(:,func_num)));    
    fprintf(' -------------------------------------------------- \n');
 end

xlabel('iteration');
ylabel('everyfit');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'-p','color','g','MarkerFaceColor','g','MarkerSize',3,'LineWidth', 0.5);  
legend('DPSO-PI','EOPSO');

% for i=1:29
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);i,f(i)
% end