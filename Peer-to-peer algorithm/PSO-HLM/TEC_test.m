%test_framaxgen
 clear all
 clc

%  mex cec17_func.cpp -DWINDOWS
global orthm best_f best_keep initial_flag fbias gbias norm_flag  %%%-----------
global accuracy                                               %%%--------------用于记录不同测试函数的最优值 / 精度有关

rand('state',sum(100*clock));%%% 伪随机树种子一样,该随机数的相似度与时间成正比，也就是说只要运算时间只要够快，随机数就相同 ，建议为rand('state',sum(100*clock)*rand(1))
warning off
% fhd=str2func('TEC_test_function');%%%函数句柄 fun,VRmin,VRmax,gbias,norm_flag,shift_flag
 fhd=str2func('cec17_func');
% fhd=str2func('benchmark_func');
fitcount=0;   % 评价次数初始化为0
D=30;        % 变量维数
N=50;% 10+floor(2*sqrt(D));        % 粒子个数
% 'TAPSO-30D'

Max_FES=1e4*D;%150000;%  % 最大评价次数。注意调整该值是要根据N的值调整maxgen的值，因为惯性权重对maxgen很敏感
shift_flag=0;
maxgen=ceil(Max_FES/N);  % 最大进化代数

norm_flag=0;  %%%-------------
runtimes=1;   % 独立运行次数
fbias=0;     %%%---------------

orthm=diag(ones(1,D));  %%%---------------
VRmax=100;
VRmin=-100;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 每个测试函数的测试范围
% VRmin=[-100,  -30,  -32,  -600,   -5.12,...    % 第二个值： CEC2008中测试函数Rosenbrock的搜索空间为 -100~100.原来的为 -32~32
%        -5.12, -500, -0.5, -2.048, -100,...
%        -100,    -pi,   -5,   -5,     -5,...
%        -10,   -100, -10,  -100,   -100,...
%        -2,    -5,   -100, -100, -100,...       
%        -50,   -100,  0,    -6.4,    12, ...
%        ...  % 下面为CEC2010测试函数的信息
%         -100,-5,-32,-100,-5,...
%         -32,-100,-100,-100,-5,...
%         -32,-100,-100,-100,-5,...
%         -32,-100,-100,-100,-100,...
%         ...  % 下面为CEC2005测试函数的信息
%         -100,-100,-100,-100,-100,...
%         -100, -600,-32,-5,-5,...
%         -0.5,-pi,-3,-100,-5,...
%         -5,-5,-5,-5,-5,...
%         -5,-5,-5,-5,-5,...
%         -6.4,-5,-5,-100,-100,...
%         -5,-5,-5,-5,-5,...
%         -6.4,-5,-100,-100];
%     
% VRmax=-VRmin;
% % VRmin(57) = 0;   % ----  CEC2005 中对这些函数进行测试时，种群初始化范围并非整个搜索空间
% VRmax(29) = 6.35;  % CEC2011-1: FM
% VRmax(30) = 60;    % Gear Train
% 
% VRmax(63) = 1;     % CEC2005 F13
% 
% VRmin(75) = 2;
% % 此3个函数的最大值为 CEC2005 中测试函数的最大值
% VRmax(63)=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
if norm_flag==1;
    VRminn=zeros(1,D);
    VRmaxn=ones(1,D);
else
    VRminn=VRmin;VRmaxn=VRmax;
end

% 若求算法求解精度时则注释掉上面的语句，采用下面的语句
accuracy=[ 1e-6, 1e2, 1e-6, 1e2, 1e-6,...
           1e-6, 1e-6, 1e-6, 0, 0,...
           0, 1e2, 0, 0, 0,...
           1e-6, 0, 0, 1e-6, -1,...
           3, -1.0316285, 1e2, 1e2, 0,...
           1e-6, 0, 0, 0, 0,...
           0, 0, 0, 0, 0,...
           0, 0, 0, 0, 0,...
           0, 0, 0, 0, 0,...
           0, 0, 0, 0, 0,...
           1e-8, 1e-8, 1e-8, 1e-8, 1e-8,...       % fun: 51-75    
           1e-8, 1e-8, 1e-8, 1e-8, 1e-8,...
           1e-8, 1e-8, 1e-8, 1e-8, 1e-8,...
           1e-8, 1e-8, 1e-8, 1e-8, 1e-8,...
           1e-8, 1e-8, 1e-8, 1e-8, 1e-8,...
           0,1e-6,1e2,1e2,0,...
           0,0,0,0,0,...
           0,0,0,0,0,...
           0,0,0,0,0,...
           0,0,0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% funchoose = [51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];  %%%要运行的测试函数的序号--函数编号
funchoose = [51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];  %%%要运行的测试函数的序号--函数编号


%   %%%%%%%%%%%%%%%%%%%%%  创建函数名数组，输出结果时用
name1='Sphere'; name2='Rosenbrock'; name3='Ackley';name4='Griewank'; name5='Rastrigin'; 
name6='Rastrigin_noncont'; name7='Schewfel'; name8='Weierstrass'; name9='EF8F2'; name10='ScafferF6';
name11='Sum of different power'; name12='Schwefel_213'; name13='Com_func1';name14='Hybrid_func1'; name15='Hybrid_func2'; 
name16='Schwefel P2.22'; name17='Step'; name18='Alpine'; name19='Salomon'; name20='Easoms';
name21='Goldstein Price'; name22='Six-hump camel back';name23='Pathologi cal';name24='Schewfel_1.2';name25='Quadric_Noise';
name26='Penalized'; name27='Schwefel_221'; name28='Fastfractal_doubledip'; name29='CEC11-01-FM'; name30='Gear_Train';
%% name31-name50 为CEC2010中的20个测试函数
name31='CEC10-F1'; name32='CEC10-F2'; name33='CEC10-F3';name34='CEC10-F4';name35='CEC10-F5';
name36='CEC10-F6'; name37='CEC10-F7'; name38='CEC10-F8';name39='CEC10-F9';name40='CEC10-F10';
name41='CEC10-F11'; name42='CEC10-F12'; name43='CEC10-F13';name44='CEC10-F14';name45='CEC10-F15';
name46='CEC10-F16'; name47='CEC10-F17'; name48='CEC10-F18';name49='CEC10-F19';name50='CEC10-F20';

%% name51-name75 为CEC2005中的25个测试函数
name51='CEC05-F1'; name52='CEC05-F2'; name53='CEC05-F3';name54='CEC05-F4';name55='CEC05-F5';
name56='CEC05-F6'; name57='CEC05-F7'; name58='CEC05-F8';name59='CEC05-F9';name60='CEC05-F10';
name61='CEC05-F11'; name62='CEC05-F12'; name63='CEC05-F13';name64='CEC05-F14';name65='CEC05-F15';
name66='CEC05-F16'; name67='CEC05-F17'; name68='CEC05-F18';name69='CEC05-F19';name70='CEC05-F20';
name71='CEC05-F21'; name72='CEC05-F22'; name73='CEC05-F23';name74='CEC05-F24';name75='CEC05-F25';

name76='frequency_modulated';name77='Shifted_nonc_rastrigin_func';name78='Shifted_Rotated_nonc_rastrigin_func';
name79='Shifted_Rotated_salomon_func'; name80='Shifted_Rotated_sphere_func';
name81='sprd_spectrum_rad_pphase'; name82='lennard_jones_potential_problem';
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
% 注意：要验证fun25，需要将shift_flag设置为1，仍选择fun24

funname=char(name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,...
            name11,name12,name13,name14,name15,name16,name17,name18,name19,name20,...
            name21,name22,name23,name24,name25,name26,name27,name28,name29,name30,...
            name31,name32,name33,name34,name35,name36,name37,name38,name39,name40,...
            name41,name42,name43,name44,name45,name46,name47,name48,name49,name50,...
            name51,name52,name53,name54,name55,name56,name57,name58,name59,name60,...
            name61,name62,name63,name64,name65,name66,name67,name68,name69,name70,...
            name71,name72,name73,name74,name75,name76,name77,name78,name79,name80,...
            name81,name82);

% 初始化保存结果的几个矩阵  %%%%%%%%%%%-----------------
CLPSO_new_fitcount_res=zeros(length(funchoose),runtimes);        % 初始化为 0      %%% 用于记录保存 每次 每个 测试函数 的 评价次数
CLPSO_new_gbestval_res=ones(length(funchoose),runtimes)*inf;     % 初始化为无穷    %%% 用于记录保存 每次 每个 测试函数 全局最优 的 值
CLPSO_new_gbest_res=ones(length(funchoose),runtimes,D)*inf;      % 初始化为无穷    %%% 用于记录保存 每次 每个 测试函数 全局最优 的 位置
bestresults=ones(runtimes,length(funchoose))*inf;   %%%  用于记录保存 每次 每个 测试函数 每次运行后 的 最优值

gbestval_res_1=ones(length(funchoose),runtimes)*inf;  %%%  用于记录保存 每次 每个 测试函数 每次运行后 的 最优值---------
bestresults_1=ones(runtimes,length(funchoose))*inf;   %%%  用于记录保存 每次 每个 测试函数 每次运行后 的 最优值--------
 
gbestval_res_2=ones(length(funchoose),runtimes)*inf;  %%%  用于记录保存 每次 每个 测试函数 每次运行后 的 最优值 --------
bestresults_2=ones(runtimes,length(funchoose))*inf; %%%  用于记录保存 每次 每个 测试函数 每次运行后 的 最优值--------

% for funnum=1:length(funchoose)  
for funnum=  	 1:1
    fun=funchoose(funnum);      %%% 对需要求解的函数 逐一进行求解，对每个函数求解 调用一下TAPSO
    
%     if fun == 76      % frequency_modulated  %%%对相应 测试函数 的 定义域 和 维度信息
%        % Max_FES=150000;
%         D=6; maxgen=ceil(Max_FES/N);  % 最大进化代数
%         VRmin(fun)=-6.4; VRmax(fun)=6.35;
%     elseif fun == 30  % Gear_Train       %%%对相应 测试函数 的 定义域 和 维度信息
%       %  Max_FES=150000;
%         D=4; maxgen=ceil(Max_FES/N);  % 最大进化代数
%         VRmin(fun)=12; VRmax(fun)=60;
%     elseif fun == 81  % sprd_spectrum_rad_pphase        %%%对相应 测试函数 的 定义域 和 维度信息
%        % Max_FES=150000;
%         D=20; maxgen=ceil(Max_FES/N);
%         VRmin(fun)=0; VRmax(fun)=2*pi;
%     elseif fun == 82  % lennard_jones_potential_problem       %%%对相应 测试函数 的 定义域 和 维度信息
%       %  Max_FES=150000;
%         D=30; maxgen=ceil(Max_FES/N);
%         
%     end

    jingdu = 0;%accuracy(fun); % 对应函数的设定可接受的精度
    initial_flag=0;
    
    timeusage=0;    %%%-------------
    fesusage=0; %%%------------
    count=0;    %%%-------------
    
    suc_times = 0;  %  满足精度要求的次数 %%%----
    for test=1:runtimes   % 每个函数独立运行 runtimes 次
        
       
%         orthm=orthm_generator(D);   % orthm为正交矩阵，目的是使seperable函数变成nonseperable函数！！！！！！！！！！！！
%         if shift_flag==1  %%%-----
% %             gbias=0.8.*(VRmin(fun)+(VRmax(fun)-VRmin(fun)).*rand(1,D));           
%             gbias=VRmin(fun)+(VRmax(fun)-VRmin(fun).*rand(1,D));  %%% gbais----------
%             % 下面几个函数是CEC2008中用到的，根据需要做了一下设置。但不影响原来的使用。
% %             if fun==2
% %                 gbias=-90+180*rand(1,D);%-1+(1+1).*rand(1,D);
% %             end
%             if fun==3
%                 gbias=-30+60*rand(1,D);
%             end
%             if fun==7
%                 gbias=-50+(0+50).*rand(1,D); % 500改为50
%             end
%             
%             if fun==1
%                 fbias=-450;
%             end
%             if fun==2
%                 fbias=390;
%             end
%             if fun==3
%                 fbias=-140;
%             end
%             if fun==4
%                 fbias=-180;
%             end
%             if fun==5
%                 fbias=-330;
%             end
%             if fun==27
%                fbias=-450; 
%             end
%         else
%             gbias=zeros(1,D);
%             fbias=0;
%         end
        
        tic; % toc  tic   计算之间程序运行所需要的时间 
         
        group_num=10;        %%%  这里要保证 group_num*group_size = ps
        group_size=3;
 
%         N=20;  % TRAPSO


        

jg=[];%++++++++++++++++++++++++++++++++++++++++++++++++
   
%         [jg,CLPSO_new_gbest,CLPSO_new_gbestval,CLPSO_new_fitcount,suc, suc_fes]= GTCAPSO_func(jingdu,fhd,maxgen,Max_FES,N,D,VRminn(fun),VRmaxn(fun),funnum,VRmin(fun),VRmax(fun),gbias,norm_flag,shift_flag); 
         [jg,CLPSO_new_gbest,CLPSO_new_gbestval,CLPSO_new_fitcount,suc, suc_fes,diversity,everyfit]= GTCAPSO_func(jingdu,fhd,maxgen,Max_FES,N,D,VRmin,VRmax,funnum); 
        t(test,funnum)=toc;                                                     %GTCAPSO_func(jingdu,fhd,Max_Gen,Max_FES,   Particle_Number,Dimension,  VRmin,      VRmax,      varargin)

       jgg(test,:)=jg(:);%++++++++++++++++++++++++++++++++++++++++++++++++
        fprintf('第 %d 次运行的最优结果为：%1.4e\n',test,CLPSO_new_gbestval(1));
        CLPSO_new_fitcount_res(funnum,test)=CLPSO_new_fitcount;   % 记录函数 fun 第 test 次运行时的总评价次数
        CLPSO_new_gbestval_res(funnum,test)=CLPSO_new_gbestval(1);   % 记录函数 fun 第 test 次运行时得到的最优解对应的结果
%         CLPSO_new_gbest_res(funnum,test,:)=CLPSO_new_gbest(1,:);       % 记录函数 fun 第 test 次运行时得到的最优解
        bestresults(test,funnum) = CLPSO_new_gbestval_res(funnum,test);
%         fitcount_used(test,fff)=CLPSO_new_fitcount_res(funnum,test);
        fprintf(' ---------------------------- \n');
        suc_times = suc_times + suc;
        if suc == 1  % 当达到设定精度时才统计其耗时           
            fesusage = fesusage + suc_fes;   % 达到精度时的 fes            
        end
        
    end  
    jggh((runtimes*funnum)+1:runtimes*(funnum+1),:)=jgg(1:runtimes,:);%++++++++++++++++++++++++++++++++++++++++++++++++
    jgg=[];%++++++++++++++++++++++++++++++++++++++++++++++++
    % 下面为计算和输出在设定求解精度条件下的成功率、所需评价次数以及运行时间
    SR(1,funnum) = suc_times/runtimes;
    if suc_times>0
        FEs(1,funnum) = fesusage/suc_times;  % 满足精度的多次运行所消耗的平均 fes。未考虑不满足精度的运行
        SP (1,funnum) = Max_FES*(1-SR(1,funnum))/SR(1,funnum) + FEs(1,funnum); % 综合评价了算法的性能：既考虑成功的，也考虑未成功的
        tu(1,funnum) = timeusage/suc_times;
    else
        FEs(1,funnum) = -1;  % 满足精度的多次运行所消耗的平均 fes。未考虑不满足精度的运行
        SP (1,funnum) = -1; % 综合评价了算法的性能：既考虑成功的，也考虑未成功的
        tu(1,funnum) = -1;
    end
    f_mean(funnum*2,:) = mean(CLPSO_new_gbestval_res(funnum,:));
     f_mean(funnum*2+1,:)=std(CLPSO_new_gbestval_res(funnum,:));
    fprintf('CEC2017\n');  % 此处根据前面调用不同的函数进行相应的输出信息
   % for i=1:length(funchoose)
        fprintf('Function %s F %d :\nAvg. fitness = %1.2e(%1.2e)\nAvg. FEs = %1.4e\n ', ...
            funname(funchoose(funnum),:),funnum,mean(CLPSO_new_gbestval_res(funnum,:)), std(CLPSO_new_gbestval_res(funnum,:)), mean(CLPSO_new_fitcount_res(funnum,:)));

        if count>0
            fprintf(' Time_Usage = %1.3e\n', timeusage/count);    % 输出达到设定精度的平均耗时
        end
        fprintf(' -------------------------------------------------- \n');
    %end

    fprintf('Avg. Suscess Ration = %1.2f\n',SR(1,funnum));
    fprintf('Avg. FEs of Success runs = %d\n', FEs(1,funnum));
    fprintf('Success Performance (SP) = %d\n', SP(1,funnum));
    fprintf('Avg. Timeusage = %1.2e\n\n', timeusage/count);% 输出达到设定精度的平均耗时
    fprintf(' =============================================================== \n\n');
    
end    
    
    
xlabel('iteration');
ylabel('diversity');
set(gca, 'Fontname', 'Times New Roman','FontSize',9);
x=1:round(6000/100):6000;
x(101)=6000;
hold on;
plot(x,log10(everyfit(x)),'-s','color','k','MarkerFaceColor','k','MarkerSize',3,'LineWidth', 0.5);  
legend('DPSO-PI','EOPSO','TAPSO','EAPSO','PSOHLM'); 

    
% end


% fprintf('Avg. time = %1.2e(%1.2e)\n', mean(tot_time1), std(tot_time1));
% mean(CLPSO_new_gbestval_res')
% CLPSO_new_gbest_res
