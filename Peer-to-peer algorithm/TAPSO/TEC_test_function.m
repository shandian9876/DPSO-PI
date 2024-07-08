function f=TEC_test_function(x,fun,VRmin,VRmax,gbias,norm_flag,shift_flag)
global orthm best_f best_keep initial_flag fbias

global jrandflag jrand lb ub  % new
persistent fhd   
persistent f_bias   % CEC2005 使用。若不测试 CEC2005 则可注释掉！！！！！！！！

[ps,D]=size(x);

if norm_flag==1
x=VRmin(1,:)+(VRmax(1,:)-VRmin(1,:)).*x;
end

if shift_flag==1
x=x-repmat(gbias,ps,1);
end

% 下面的应该是每个函数的理论最优值
greal=[0, 1, 0, 0, 0,...                       % fun: 1-30
       0, 4.209687462275036e+002, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, -1,...
       3, -1.031628,0,0, 0,...
       0,-1,-1,-1,-1,...
       0, 0, 0, 0, 0,...       % fun: 31-50
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       0, 0, 0, 0, 0,...
       -450, -450, -450, -450, -310,...       % fun: 51-75    
       390, -180, -140, -330, -330,...
       90, -460, -130, -300, 120,...
       120, 120, 10, 10, 10,...
       360, 360, 360, 260, 260,...
       0,-330,-330,0,0];
   
% x=x-greal(fun);  % 此句与line39句一起似乎没啥意义啊？
% x=x*orthm;
% x=x+greal(fun);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 下面为常用的一些函数
if fun==1
    %sphere with noise
%     f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
%     f=sum(x.^2,2);
    f=sum(x.^2);%+fbias;

elseif fun==2
    %rosenbrock
    f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);%+fbias;
    
    
elseif fun==3
    %ackley
    f=sum(x.^2,2);
    f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);%+fbias;
    
elseif fun==4
    %griewank
    f=1;
    for i=1:D
        f=f.*cos(x(:,i)./sqrt(i));
    end
    f=sum(x.^2,2)./4000-f+1;%+fbias;

elseif fun==5
    %rastrigin
    f=sum(x.^2-10.*cos(2.*pi.*x)+10,2)+fbias;

elseif fun==6
    %rastrigin_noncont
    x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
%     x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2));
    f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);%+fbias;
    
elseif fun==7
    %schewfel
    f = 418.9828872724338*D - sum(x.*sin(sqrt(abs(x))),2);%+fbias;
    
elseif fun==8
    %weierstrass
    x=x+0.5;
    a = 0.5;
    b = 3;
    kmax = 20;
    c1(1:kmax+1) = a.^(0:kmax);
    c2(1:kmax+1) = 2*pi*b.^(0:kmax);
    f=0;
    c=-w(0.5,c1,c2);
    for i=1:D
    f=f+w(x(:,i)',c1,c2);
    end
    f=f+c*D;%+fbias;
    
elseif fun==9
    %EF8F2
    f=0;
    for i=1:(D-1)
        f=f+F8F2(x(:,[i,i+1]));
    end
    f=f+F8F2(x(:,[D,1]));%+fbias;
    
elseif fun==10
    %E_ScafferF6
    fhd=str2func('ScafferF6');
    f=0;
    for i=1:(D-1)
        f=f+feval(fhd,(x(:,i:i+1)));
    end
    f=f+feval(fhd,x(:,[D,1]));%+fbias;
    
elseif fun==11  % Sum of different power
    f=0;
    for i=1:D
        f=f+abs(x(:,i)).^(i+1);
    end
    f=f;%+fbias;

elseif fun==12
    f=schwefel_213(x);%+fbias;

elseif fun==13
    f=com_func1(x);%+fbias;   

elseif fun==14
    f=hybrid_func1(x);%+fbias;
    
elseif fun==15
    f=hybrid_func2(x);%+fbias;
    
elseif fun==16   % 自己加的： Schwefel P2.22
    f1=sum(abs(x),2);
    f2=prod(abs(x),2);
    f=f1+f2;%+fbias;
    
elseif fun==17   % 自己加的： Step
    f=sum((floor(x+0.5)).^2,2);%+fbias;

elseif fun==18   % 自己加的：Alpine
    f=sum(abs(x.*sin(x)+0.1.*x),2);%+fbias; 
    
elseif fun==19   % 自己加的：Salomon
    persistent o M
    [ps,D]=size(x);
    if initial_flag==0
     load sphere_func_data   %  CEC2005中，只支持10D，30D,50D
%    load sphere_shift_func_data   %  CEC2008中，只支持100D，500D,1000D
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;    
    c=1;   %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
        if D==2,load elliptic_M_D2,  
        elseif D==10,load elliptic_M_D10,  
        elseif D==30,load elliptic_M_D30, 
        elseif D==50,load elliptic_M_D50,  
        elseif D==100,load elliptic_M_D100,
        else 
            M=rot_matrix(D,c);
            M=M.*(1+0.3.*normrnd(0,1,D,D));
        end
        %----------------------------------------------------
    end
    x=x-repmat(o,ps,1);
%     x=x*M; %---------------------------------------------------
    f1=2*pi.*sqrt(sum(x.^2,2));
    f2=0.1.*sqrt(sum(x.^2,2))+1;
    f=-cos(f1)+f2;%+fbias;
    
elseif fun==20   % 自己加的：Easom's
    f1=(x(1)-pi)^2;
    f2=(x(2)-pi)^2;
    f=-cos(x(1))*cos(x(2))*exp(-(f1+f2));%+fbias;
    
elseif fun==21   % 自己加的：Goldstein Price's
    f=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2))* (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));+fbias;

elseif fun==22   % 自己加的：Six-hump camel back
    f=(4-2.1*x(1)^2+x(1)^4/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;%+fbias;
    
elseif fun==23   % 自己加的：Pathological
    f=0;    
    for j=1:D-1
        y1=(sin(sqrt((x(:,j+1)*x(:,j+1))+100*x(:,j)*x(:,j))))^2-0.5;
        y2=0.001*((x(:,j+1)*x(:,j+1)-2*x(:,j+1)*x(:,j)+x(:,j)*x(:,j)))^2+1;
        f=f+y1/y2+0.5;
    end
    f=f;%+fbias;

elseif fun==24   % 自己加的：Quadric    
    f = 0;  
    xj=0;
    for j=1:D
        xj = xj + sum(x(:,1:j),2).^2;
    end
    f = xj;
%  
%     f=0; 
%     y2=0;   
%     for i=1:D 
%         y1=0;
%         for j=1:i
%             y1=y1+x(:,j);
%         end
%         y2=y2+y1*y1;
%     end
%     f=y2+fbias;
    
elseif fun==25   % 自己加的：Quadric_Noise    
    f=0; 
    y2=0;  
    x=x+rand(1,D);
    for i=1:D 
        y1=0;
        for j=1:i
            y1=y1+x(:,j);
        end
        y2=y2+y1*y1;
    end
    f=y2;%+fbias;
    
elseif fun==26   % 自己加的：Penalized
    f=0;
    xu=u(x,10,100,4);
    y=1+(x+1)/4;
    y1=10*(sin(pi*y(1))*sin(pi*y(1)))+(y(D)-1)*(y(D)-1) ;
    y2=0;
    for kk=1:D-1
       y2=y2+(y(kk)-1)*(y(kk)-1)*(1+10*sin(pi*y(kk+1))*sin(pi*y(kk+1))) ;
    end
    y3=sum(xu);
    f=pi/D*(y1+y2)+y3;%+fbias;
    
elseif fun==27   % 自己加的：Schwefel's P 2.21
    f=0;
    f = max(abs(x), [], 2);
    
elseif fun==28   % 自己加的：fastfractal_doubledip
%     f=0;
%     f = max(abs(x), [], 2);

elseif fun==29   % 自己加的：CEC11-01
    theta=2*pi/100;
    f=0;
   for t=0:100
       y_t=x(1)*sin(x(2)*t*theta+x(3)*sin(x(4)*t*theta+x(5)*sin(x(6)*t*theta)));
       y_0_t=1*sin(5*t*theta-1.5*sin(4.8*t*theta+2*sin(4.9*t*theta)));
       f=f+(y_t-y_0_t)^2;
   end
   
elseif fun==30   % 自己加的：Gear_Train    
   f=(1/6.931-(x(1)*x(2))/(x(3)*x(4)))*(1/6.931-(x(1)*x(2))/(x(3)*x(4)));    
end


%%  fun31-fun50 为CEC2010的测试函数。当不使用这些函数时，请将line41-line840注释掉！！！
% % Separable, D = 1000
% %   my_m 和 dim 主要用于 FUN4~FUN18. 因为这几个函数原版本只支持1000维，这里做了一些修改
% %   此外，后面每一函数中的语句：o = o(1:D); 也系自己添加，目的是将1000维的 o 调整为自己需要的 o  .
%     
% 	if (fun ==  31) fhd = str2func('elliptic_shift_func');
%     elseif (fun ==  32) fhd = str2func('rastrigin_shift_func');
%     elseif (fun ==  33) fhd = str2func('ackley_shift_func');
%     % Single-group m-nonseparable, D = 1000, m = 50
% 	elseif (fun ==  34) fhd = str2func('elliptic_group1_shift_rot_func');
%     elseif (fun ==  35) fhd = str2func('rastrigin_group1_shift_rot_func');
%     elseif (fun ==  36) fhd = str2func('ackley_group1_shift_rot_func');
%     elseif (fun ==  37) fhd = str2func('schwefel_group1_shift_func');
%     elseif (fun ==  38) fhd = str2func('rosenbrock_group1_shift_func');
%     % D/(2m)-group m-nonseparable, D = 1000, m = 50
% 	elseif (fun ==  39) fhd = str2func('elliptic_group10_shift_rot_func');
%     elseif (fun == 40) fhd = str2func('rastrigin_group10_shift_rot_func');
%     elseif (fun == 41) fhd = str2func('ackley_group10_shift_rot_func');
%     elseif (fun == 42) fhd = str2func('schwefel_group10_shift_func');
%     elseif (fun == 43) fhd = str2func('rosenbrock_group10_shift_func');
%     % D/m-group m-nonseparable, D = 1000, m = 50
% 	elseif (fun == 44) fhd = str2func('elliptic_group20_shift_rot_func');
%     elseif (fun == 45) fhd = str2func('rastrigin_group20_shift_rot_func');
%     elseif (fun == 46) fhd = str2func('ackley_group20_shift_rot_func');
%     elseif (fun == 47) fhd = str2func('schwefel_group20_shift_func');
%     elseif (fun == 48) fhd = str2func('rosenbrock_group20_shift_func');
%     % Fully-nonseparable, D = 1000
% 	elseif (fun == 49) fhd = str2func('schwefel_shift_func');
%     elseif (fun == 50) fhd = str2func('rosenbrock_shift_func');
%     end
%     
%     if fun>=31 && fun<=50
%         jrandflag = 0;  
%         if (jrandflag == 1)    
%             jrand = Randomizer(func_num); 
%         end
%     end
%     
%     if fun>=31 && fun<=50
%        f = feval(fhd, x); 
%     end
%     
%     
%     
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sphere Function 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = sphere_func(x)
% 
% fit = sum(x.*x, 2);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Elliptic Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_func(x)
% 
% [ps, D] = size(x);
% a = 1e+6;
% fit = 0;
% for i = 1:D
%    fit = fit+a.^((i-1)/(D-1)).*x(:,i).^2;
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rotated Elliptic Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_rot_func(x, M)
% 
% x = x*M;
% fit = elliptic_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Schwefel's Problem 1.2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = schwefel_func(x)
% 
% [ps D] = size(x);
% fit = 0;
% for i = 1:D
% 	fit = fit + sum(x(:,1:i),2).^2;
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rosenbrock's Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rosenbrock_func(x)
% 
% [ps D] = size(x);
% fit = sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rastrigin's Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_func(x)
% 
% fit = sum(x.*x-10*cos(2*pi*x)+10, 2);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rotated Rastrigin's Fucntion 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_rot_func(x, M)
% 
% x = x*M;
% fit = rastrigin_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Ackley's Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_func(x)
% 
% [ps, D] = size(x);
% fit = sum(x.^2,2);
% fit = 20-20.*exp(-0.2.*sqrt(fit./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rotated Ackley's Function 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_rot_func(x, M)
% 
% x = x*M;
% fit = ackley_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F1: Shifted Elliptic Function
% % D = 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o
% 
% [ps D] = size(x);
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
% 		save 'datafiles/f01_o.mat' o;
%     else
% 		load 'datafiles/f01_o.mat';
% 		o = o(1:D);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = elliptic_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F2: Shifted Rastrigin's Function
% % D = 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o
% 
% [ps D] = size(x);
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
% 		save 'datafiles/f02_o.mat' o;
%     else
% 		load 'datafiles/f02_o.mat';
% 		o = o(1:D);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = rastrigin_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F3: Shifted Ackley's Function
% % D = 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o
% 
% [ps D] = size(x);
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
% 		save 'datafiles/f03_o.mat' o;
%     else
% 		load 'datafiles/f03_o.mat' o;
% 		o = o(1:D);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = ackley_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F4: Single Group Shifted and Rotated Elliptic Function
% % D = 1000, m = 50
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_group1_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% %  my_m = 5;    % 此处为自己修改，CEC2010中设置的是 m=50
% %  dim = 30;    % 此处为变量维数30，CEC2010中设置的变量维数为 1000 
%  
% [ps D] = size(x);
% % m = my_m;
% c = 1;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
%         M = jrand.createRotMatrix(m);
% 		save 'datafiles/f04_opm.mat' o p M;
%     else
% 		load 'datafiles/f04_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F4 error: only support D = 1000 now');
%         exit(4);
%     end
% 	initial_flag = 1;
%     o = o(1:D);
%     p = randperm(D);
%     M = rot_matrix(m,c);
% end
% 
% a = 1e+6;
% x = x-repmat(o,ps,1);
% fit = a*elliptic_rot_func(x(:,p(1:m)), M) + elliptic_func(x(:,p((m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F5: Single Group Shifted and Rotated Rastrigin's Function
% % D = 1000, m = 50
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_group1_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = 50;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f05_opm.mat' o p M;
%     else
% 		load 'datafiles/f05_opm.mat';
%     end
%     if (D ~= 1000)%(D ~= dim)
%         disp('F5 error: only support D = 1000 now');
%         exit(5);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% a = 1e+6;
% x = x-repmat(o,ps,1);
% fit = a*rastrigin_rot_func(x(:,p(1:m)), M) + rastrigin_func(x(:,p((m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F6: Single Group Shifted and Rotated Ackley's Function
% % D = 1000, m = 50
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_group1_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f06_opm.mat' o p M;
%     else
% 		load 'datafiles/f06_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F6 error: only support D = 1000 now');
%         exit(6);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% a = 1e+6;
% x = x-repmat(o,ps,1);
% fit = a*ackley_rot_func(x(:,p(1:m)), M) + ackley_func(x(:,p((m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F7: Single Group Shifted Schwefel's Problem 1.2
% % D = 1000, m = 50
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = schwefel_group1_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = 50;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f07_op.mat' o p;
%     else
% 		load 'datafiles/f07_op.mat';
%     end
%     if (D ~= 1000)
%         disp('F7 error: only support D = 1000 now');
%         exit(7);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% a = 1e+6;
% x = x-repmat(o,ps,1);
% fit = a*schwefel_func(x(:,p(1:m))) + sphere_func(x(:,p((m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F8: Single Group Shifted Rosenbrock's Function
% % D = 1000, m = 50
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rosenbrock_group1_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = my_m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub-1);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f08_op.mat' o p;
%     else
% 		load 'datafiles/f08_op.mat';
%     end
%     if (D ~= dim)
%         disp('F8 error: only support D = 1000 now');
%         exit(8);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% a = 1e+6;
% x = x-repmat(o,ps,1);
% fit = a*rosenbrock_func(x(:,p(1:m))) + sphere_func(x(:,p((m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F9: D/(2m)-group Shifted and Rotated Elliptic Function
% % D = 1000, m = 50, D/(2m) = 10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_group10_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m/2;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f09_opm.mat' o p M;
%     else
% 		load 'datafiles/f09_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F9 error: only support D = 1000 now');
%         exit(9);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + elliptic_rot_func(x(:,p(index)), M);
% end
% fit = fit + elliptic_func(x(:,p((G*m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F10: D/(2m)-group Shifted and Rotated Rastrigin's Function
% % D = 1000, m = 50, D/(2m) = 10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_group10_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m/2;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f10_opm.mat' o p M;
%     else
% 		load 'datafiles/f10_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F10 error: only support D = 1000 now');
%         exit(10);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + rastrigin_rot_func(x(:,p(index)), M);
% end
% fit = fit + rastrigin_func(x(:,p((G*m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F11: D/(2m)-group Shifted and Rotated Ackley's Function
% % D = 1000, m = 50, D/(2m) = 10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_group10_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m/2;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f11_opm.mat' o p M;
%     else
% 		load 'datafiles/f11_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F11 error: only support D = 1000 now');
%         exit(11);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + ackley_rot_func(x(:,p(index)), M);
% end
% fit = fit + ackley_func(x(:,p((G*m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F12: D/(2m)-group Shifted Schwefel's Problem 1.2
% % D = 1000, m = 50, D/(2m) = 10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = schwefel_group10_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m/2;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f12_op.mat' o p;
%     else
% 		load 'datafiles/f12_op.mat';
%     end
%     if (D ~= dim)
%         disp('F12 error: only support D = 1000 now');
%         exit(12);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + schwefel_func(x(:,p(index)));
% end
% fit = fit + sphere_func(x(:,p((G*m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F13: D/(2m)-group Shifted Rosenbrock's Function
% % D = 1000, m = 50, D/(2m) = 10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rosenbrock_group10_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m/2;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub-1);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f13_op.mat' o p;
%     else
% 		load 'datafiles/f13_op.mat';
%     end
%     if (D ~= dim)
%         disp('F13 error: only support D = 1000 now');
%         exit(13);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + rosenbrock_func(x(:,p(index)));
% end
% fit = fit + sphere_func(x(:,p((G*m+1):end)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F14: D/m-group Shifted and Rotated Elliptic Function
% % D = 1000, m = 50, D/m = 20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = elliptic_group20_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f14_opm.mat' o p M;
%     else
% 		load 'datafiles/f14_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F14 error: only support D = 1000 now');
%         exit(14);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + elliptic_rot_func(x(:,p(index)), M);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F15: D/m-group Shifted and Rotated Rastrigin's Function
% % D = 1000, m = 50, D/m = 20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rastrigin_group20_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f15_opm.mat' o p M;
%     else
% 		load 'datafiles/f15_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F15 error: only support D = 1000 now');
%         exit(15);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + rastrigin_rot_func(x(:,p(index)), M);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F16: D/m-group Shifted and Rotated Ackley's Function
% % D = 1000, m = 50, D/m = 20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = ackley_group20_shift_rot_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p M
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		M = jrand.createRotMatrix(m);
% 		save 'datafiles/f16_opm.mat' o p M;
%     else
% 		load 'datafiles/f16_opm.mat';
%     end
%     if (D ~= dim)
%         disp('F16 error: only support D = 1000 now');
%         exit(16);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + ackley_rot_func(x(:,p(index)), M);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F17: D/m-group Shifted Schwefel's Problem 1.2
% % D = 1000, m = 50, D/m = 20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = schwefel_group20_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f17_op.mat' o p;
%     else
% 		load 'datafiles/f17_op.mat';
%     end
%     if (D ~= dim)
%         disp('F17 error: only support D = 1000 now');
%         exit(17);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + schwefel_func(x(:,p(index)));
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F18: D/m-group Shifted Rosenbrock's Function
% % D = 1000, m = 50, D/m = 20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rosenbrock_group20_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o p
% 
% [ps D] = size(x);
% m = my_m;
% G = D/m;
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub-1);
%         o = o';
%         p = jrand.createPermVector(D);
%         p = p'+1;
% 		save 'datafiles/f18_op.mat' o p;
%     else
% 		load 'datafiles/f18_op.mat';
%     end
%     if (D ~= dim)
%         disp('F18 error: only support D = 1000 now');
%         exit(18);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = 0;
% for k = 1:G
%     index = ((k-1)*m+1):(k*m);
%     fit = fit + rosenbrock_func(x(:,p(index)));
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F19: Shifted Schwefel's Problem 1.2
% % D = 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = schwefel_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o
% 
% [ps D] = size(x);
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub);
%         o = o';
% 		save 'datafiles/f19_o.mat' o;
%     else
% 		load 'datafiles/f19_o.mat';
%     end
%     if (D ~= dim)
%         disp('F19 error: only support D = 1000 now');
%         exit(19);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = schwefel_func(x);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % F20: Shifted Rosenbrock's Function
% % D = 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fit = rosenbrock_shift_func(x)
% global initial_flag jrandflag jrand lb ub
% persistent o
% 
% [ps D] = size(x);
% if (initial_flag == 0)
%     if (jrandflag == 1)
%         o = jrand.createShiftVector(D, lb, ub-1);
%         o = o';
% 		save 'datafiles/f20_o.mat' o;
%     else
% 		load 'datafiles/f20_o.mat';
%     end
%     if (D ~= dim)
%         disp('F20 error: only support D = 1000 now');
%         exit(20);
%     end
% 	initial_flag = 1;
% end
% o = o(1:D);
% x = x-repmat(o,ps,1);
% fit = rosenbrock_func(x);
% 
% 
% function f=F8F2(x)
% f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
% f=1+f2.^2./4000-cos(f2);
% 
% function f=ScafferF6(x)
% f=0.5+(sin(sqrt(x(:,1).^2+x(:,2).^2)).^2-0.5)./(1+0.001*(x(:,1).^2+x(:,2).^2)).^2;
% 
%     % 	12.Schwefel's Problem 2.13
% function f=schwefel_213(x)%after Fletcher and Powell
% global initial_flag
% persistent a b A alpha
% [ps,D]=size(x);
% if initial_flag==0
%     initial_flag=1;
%     load schwefel_213_data
%     if length(alpha)>=D
%         alpha=alpha(1:D);a=a(1:D,1:D);b=b(1:D,1:D);
%     else
%         alpha=-3+6*rand(1,D);
%         a=round(-100+200.*rand(D,D));
%         b=round(-100+200.*rand(D,D));
%     end
%     alpha=repmat(alpha,D,1);
%     A=sum(a.*sin(alpha)+b.*cos(alpha),2);
% end
% 
% for i=1:ps
%     xx=repmat(x(i,:),D,1);
%     B=sum(a.*sin(xx)+b.*cos(xx),2);
%     f(i,1)=sum((A-B).^2,1);
% end
% 
% %---------------------------------------------------
% %   1.com Composition Function 1
% function fit=com_func1(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load com_func1_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=zeros(1,D);
%     func.f1=str2func('fsphere');
%     func.f2=str2func('fsphere');
%     func.f3=str2func('fsphere');
%     func.f4=str2func('fsphere');
%     func.f5=str2func('fsphere');
%     func.f6=str2func('fsphere');
%     func.f7=str2func('fsphere');
%     func.f8=str2func('fsphere');
%     func.f9=str2func('fsphere');
%     func.f10=str2func('fsphere');
%     bias=((1:fun_num)-1).*100;
%     sigma=ones(1,fun_num);
%     lamda=5/100.*ones(fun_num,1);
%     lamda=repmat(lamda,1,D);
%     for i=1:fun_num
%         eval(['M.M' int2str(i) '=diag(ones(1,D));']);
%     end
% end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
% 
% %---------------------------------------------------
% %   2.com Composition Function 2
% function fit=com_func2(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load com_func2_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=zeros(1,D);
%     func.f1=str2func('fgriewank');
%     func.f2=str2func('fgriewank');
%     func.f3=str2func('fgriewank');
%     func.f4=str2func('fgriewank');
%     func.f5=str2func('fgriewank');
%     func.f6=str2func('fgriewank');
%     func.f7=str2func('fgriewank');
%     func.f8=str2func('fgriewank');
%     func.f9=str2func('fgriewank');
%     func.f10=str2func('fgriewank');
%     bias=((1:fun_num)-1).*100;
%     sigma=ones(1,fun_num);
%     lamda=5/100.*ones(fun_num,1);
%     lamda=repmat(lamda,1,D);
%     if D==10
%     load com_func2_M_D10,
%     elseif D==30
%     load com_func2_M_D30,
%     else
%         for i=1:fun_num
%             eval(['M.M' int2str(i) '=orthm_generator(D);']);
%         end
%     end
% end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%     
% %---------------------------------------------------
% %   3.com Composition Function 3
% function fit=com_func3(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load com_func3_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=zeros(1,D);
%     func.f1=str2func('frastrigin');
%     func.f2=str2func('frastrigin');
%     func.f3=str2func('frastrigin');
%     func.f4=str2func('frastrigin');
%     func.f5=str2func('frastrigin');
%     func.f6=str2func('frastrigin');
%     func.f7=str2func('frastrigin');
%     func.f8=str2func('frastrigin');
%     func.f9=str2func('frastrigin');
%     func.f10=str2func('frastrigin');
%     bias=((1:fun_num)-1).*100;
%     sigma=ones(1,fun_num);
%     lamda=ones(fun_num,1);
%     lamda=repmat(lamda,1,D);
%     if D==10
%     load com_func3_M_D10,
%     elseif D==30
%     load com_func3_M_D30,
%     else
%         for i=1:fun_num
%             eval(['M.M' int2str(i) '=orthm_generator(D);']);
%         end
%     end
% end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
% 
% 
% %----------------------------------------------------------------
% %   4.	Rotated Hybrid Composition Function 1
% function fit=hybrid_func1(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load hybrid_func1_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=0;
%     func.f1=str2func('fackley');
%     func.f2=str2func('fackley');
%     func.f3=str2func('frastrigin');
%     func.f4=str2func('frastrigin');
%     func.f5=str2func('fweierstrass');
%     func.f6=str2func('fweierstrass');
%     func.f7=str2func('fgriewank');
%     func.f8=str2func('fgriewank');
%     func.f9=str2func('fsphere');
%     func.f10=str2func('fsphere');
%     bias=((1:fun_num)-1).*100;
%     sigma=ones(1,fun_num);
%     lamda=[5/32; 5/32; 1; 1; 10; 10; 5/100; 5/100;  5/100; 5/100];
%     lamda=repmat(lamda,1,D);
%     if D==10
%     load hybrid_func1_M_D10,
%     elseif D==30
%     load hybrid_func1_M_D30,
%     else
%         for i=1:fun_num
%             eval(['M.M' int2str(i) '=orthm_generator(D);']);
%         end
%     end
% end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
% %----------------------------------------------------------------
% %   5.Rotated Hybrid Composition Function 2	
% function fit=hybrid_func2(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load hybrid_func2_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=zeros(1,D);
%     func.f1=str2func('frastrigin');
%     func.f2=str2func('frastrigin');
%     func.f3=str2func('fweierstrass');
%     func.f4=str2func('fweierstrass');
%     func.f5=str2func('fgriewank');
%     func.f6=str2func('fgriewank');
%     func.f7=str2func('fackley');
%     func.f8=str2func('fackley');
%     func.f9=str2func('fsphere');
%     func.f10=str2func('fsphere');
%     bias=((1:fun_num)-1).*100;
%     sigma=ones(1,fun_num); 
%     lamda=[1/5;1/5;10;10;5/100;5/100;5/32;5/32;5/100;5/100];
%     lamda=repmat(lamda,1,D);
%     if D==10
%     load hybrid_func2_M_D10,
%     elseif D==30
%     load hybrid_func2_M_D30,
%     else
%         for i=1:fun_num
%             eval(['M.M' int2str(i) '=orthm_generator(D);']);
%         end
%     end
%     end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
% %---------------------------------------------------------------------
% %   6.	Rotated Hybrid Composition Function 3
% function fit=hybrid_func3(x)
% global initial_flag
% persistent  fun_num func o sigma lamda bias M
% if initial_flag==0
%     [ps,D]=size(x);
%     initial_flag=1;
%     fun_num=10;
%     load hybrid_func2_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
%          o=-5+10*rand(fun_num,D);
%     end
%     o(10,:)=0;
%     func.f1=str2func('frastrigin');
%     func.f2=str2func('frastrigin');
%     func.f3=str2func('fweierstrass');
%     func.f4=str2func('fweierstrass');
%     func.f5=str2func('fgriewank');
%     func.f6=str2func('fgriewank');
%     func.f7=str2func('fackley');
%     func.f8=str2func('fackley');
%     func.f9=str2func('fsphere');
%     func.f10=str2func('fsphere');
%     bias=((1:fun_num)-1).*100;
%     sigma=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
%     lamda=[1/5;1/5;10;10;5/100;5/100;5/32;5/32;5/100;5/100];lamda=lamda.*sigma';
%     lamda=repmat(lamda,1,D);
%     if D==10
%     load hybrid_func2_M_D10,
%     elseif D==30
%     load hybrid_func2_M_D30,
%     else
%         for i=1:fun_num
%             eval(['M.M' int2str(i) '=orthm_generator(D);']);
%         end
%     end
% end
% fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
% %----------------------------------
% function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M)
% [ps,D]=size(x);
% for i=1:fun_num
%     oo=repmat(o(i,:),ps,1);
%     weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2));
% end
% 
% [tmp,tmpid]=sort(weight,2);
% for i=1:ps
%     weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
% end
% weight=weight./repmat(sum(weight,2),1,fun_num);
% 
% fit=0;
% for i=1:fun_num
%     oo=repmat(o(i,:),ps,1);
%     eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1))*M.M' int2str(i) ');']);
%     x1=5*ones(1,D);
%     eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))*M.M' int2str(i) ');']);
%     fit1=2000.*f./f1;
%     fit=fit+weight(:,i).*(fit1+bias(i));
% end
% %-------------------------------------------------
% %basic functions
% 
% function f=fsphere(x)
% %Please notice there is no use to rotate a sphere function, with rotation
% %here just for a similar structure as other functions and easy programming
% [ps,D]=size(x);
% f=sum(x.^2,2);
% %--------------------------------
% function f=fsphere_noise(x)
% [ps,D]=size(x);
% f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
% %--------------------------------
% function f=fgriewank(x)
% [ps,D]=size(x);
% f=1;
% for i=1:D
%     f=f.*cos(x(:,i)./sqrt(i));
% end
% f=sum(x.^2,2)./4000-f+1;
% %--------------------------------
% function f=fackley(x)
% [ps,D]=size(x);
% f=sum(x.^2,2);
% f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
% %--------------------------------
% function f=frastrigin(x)
% [ps,D]=size(x);
% f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
% %--------------------------------
% function f=frastrigin_noncont(x)
% [ps,D]=size(x);
% x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
% f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
% %--------------------------------
% function [f]=fweierstrass(x)
% [ps,D]=size(x);
% x=x+0.5;
% a = 0.5;
% b = 3;
% kmax = 20;
% c1(1:kmax+1) = a.^(0:kmax);
% c2(1:kmax+1) = 2*pi*b.^(0:kmax);
% f=0;
% c=-w(0.5,c1,c2);
% for i=1:D
% f=f+w(x(:,i)',c1,c2);
% end
% f=f+c*D;
% 
% function y = w(x,c1,c2)
% y = zeros(length(x),1);
% for k = 1:length(x)
% 	y(k) = sum(c1 .* cos(c2.*x(:,k)));
% end
% %--------------------------------
% function f=fE_ScafferF6(x)
% fhd=str2func('ScafferF6');
% [ps,D]=size(x);
% 
% f=0;
% for i=1:(D-1)
%     f=f+feval(fhd,(x(:,i:i+1)));
% end
%     f=f+feval(fhd,x(:,[D,1]));
% %--------------------------------    
% function f=fE_ScafferF6_noncont(x)
% fhd=str2func('ScafferF6');
% [ps,D]=size(x);
% x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
% f=0;
% for i=1:(D-1)
%     f=f+feval(fhd,(x(:,i:i+1)));
% end
%     f=f+feval(fhd,x(:,[D,1]));
% %------------------------------
% function f=fEF8F2(x)
% [ps,D]=size(x);
% f=0;
% for i=1:(D-1)
%     f=f+F8F2(x(:,[i,i+1]));
% end
%     f=f+F8F2(x(:,[D,1]));
% 
% %--------------------------------
% function f=fschwefel_102(x)
% [ps,D]=size(x);
% f=0;
% for i=1:D
%     f=f+sum(x(:,1:i),2).^2;
% end
% %--------------------------------
% function f=felliptic(x)
% [ps,D]=size(x);
% a=1e+6;
% f=0;
% for i=1:D
% f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
% end
% %--------------------------------
% 
% % classical Gram Schmid 
%  function [q,r] = cGram_Schmidt (A)
% % computes the QR factorization of $A$ via
% % classical Gram Schmid 
% % 
%  [n,m] = size(A); 
%  q = A;    
%  for j=1:m
%      for i=1:j-1 
%          r(i,j) = q(:,j)'*q(:,i);
%      end
%      for i=1:j-1   
%        q(:,j) = q(:,j) -  r(i,j)*q(:,i);
%      end
%      t =  norm(q(:,j),2 ) ;
%      q(:,j) = q(:,j) / t ;
%      r(j,j) = t  ;
%  end 
%  
% function M=rot_matrix(D,c)
% A=normrnd(0,1,D,D);
% P=cGram_Schmidt(A);
% A=normrnd(0,1,D,D);
% Q=cGram_Schmidt(A);
% u=rand(1,D);
% D=c.^((u-min(u))./(max(u)-min(u)));
% D=diag(D);
% M=P*D*Q;
% 
% %%%%%%%% 下面的函数 ff 仅作为切断与后面“if fun==1”之间的联系，无其它意义。
% %%%%%%%% 若无 ff 函数，则在优化上面的 F20: Shifted Rosenbrock's Function 时会继续执行
% %%%%%%%% 后面的 “if fun==1”语句，从而出错。
% function fit = ff(x)   
% fit = -1;
%% ------------- 以上为 CEC2010 测试函数 --------------------------
% 

%%  fun51-fun75 为CEC2005的测试函数。-----------------------------------
func_num = fun-50;
if initial_flag==0
    if func_num==1 fhd=str2func('sphere_func'); %[-100,100]
    elseif func_num==2 fhd=str2func('schwefel_102'); %[-100,100]
    elseif func_num==3 fhd=str2func('high_cond_elliptic_rot_func'); %[-100,100]
    elseif func_num==4 fhd=str2func('schwefel_102_noise_func'); %[-100,100]
    elseif func_num==5 fhd=str2func('schwefel_206'); %[no bound],initial[-100,100];
    elseif func_num==6 fhd=str2func('rosenbrock_func'); %[-100,100]
    elseif func_num==7 fhd=str2func('griewank_rot_func'); %[-600,600]
    elseif func_num==8 fhd=str2func('ackley_rot_func'); %[-32,32]
    elseif func_num==9 fhd=str2func('rastrigin_func'); %[-5,5]
    elseif func_num==10 fhd=str2func('rastrigin_rot_func'); %[-5,5]
    elseif func_num==11 fhd=str2func('weierstrass_rot'); %[-0.5,0.5]
    elseif func_num==12 fhd=str2func('schwefel_213'); %[-pi,pi]
    elseif func_num==13 fhd=str2func('EF8F2_func'); %[-3,1] 
    elseif func_num==14 fhd=str2func('E_ScafferF6_func'); %[-100,100]
    elseif func_num==15 fhd=str2func('hybrid_func1'); %[-5,5]
    elseif func_num==16 fhd=str2func('hybrid_rot_func1'); %[-5,5]
    elseif func_num==17 fhd=str2func('hybrid_rot_func1_noise'); %[-5,5]        
    elseif func_num==18 fhd=str2func('hybrid_rot_func2'); %[-5,5]    
    elseif func_num==19 fhd=str2func('hybrid_rot_func2_narrow'); %[-5,5]    
    elseif func_num==20 fhd=str2func('hybrid_rot_func2_onbound'); %[-5,5]  
    elseif func_num==21 fhd=str2func('hybrid_rot_func3'); %[-5,5]    
    elseif func_num==22 fhd=str2func('hybrid_rot_func3_highcond'); %[-5,5]  
    elseif func_num==23 fhd=str2func('hybrid_rot_func3_noncont'); %[-5,5]   
    elseif func_num==24 fhd=str2func('hybrid_rot_func4'); %[-5,5]  
    elseif func_num==25 fhd=str2func('hybrid_rot_func4'); %[-5,5]  
    elseif func_num ==26 fhd = str2func('frequency_modulated');
        
    elseif func_num ==27 fhd = str2func('nonc_rastrigin_func');     %[-5,5] 自己加的 shifted Nonc.Rastrigin
    elseif func_num ==28 fhd = str2func('nonc_rastrigin_rot_func'); %[-5,5] 自己加的 shifted rotated Nonc.Rastrigin
    elseif func_num == 29 fhd = str2func('Shifted_Rotated_salomon_func'); % [-100,100]  自己加的 Shifted_Rotated_salomon
    elseif func_num == 30 fhd = str2func('Shifted_Rotated_sphere_func'); % [-100,100]  自己加的 Shifted_Rotated_salomon
    elseif func_num == 31 fhd = str2func('sprd_spectrum_rad_pphase'); % [0,2*pi]  自己加的 CEC2011-F7
    elseif func_num == 32 fhd = str2func('lennard_jones_potential_problem'); % 自己加的 CEC2011-F2
    end
    load fbias_data;    
end

% fn = func_num - 50;  % 因为 CEC2005 的函数被添加到本代码中是从第 51 个编号开始的，所以要：fun-50
if fun>=51 
% f=feval(fhd,x)+f_bias(func_num);
% f=f-f_bias(func_num);
f=feval(fhd,x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Unimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	1.Shifted Sphere Function 
function fit=sphere_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
     load sphere_func_data   %  CEC2005中，只支持10D，30D,50D
%    load sphere_shift_func_data   %  CEC2008中，只支持100D，500D,1000D
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
    
%     c=1;   %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
%     %---------------------------------------------------
%     if D==2,load rosenbrock_M_D2,
%     elseif D==10,load rosenbrock_M_D10,
%     elseif D==30,load rosenbrock_M_D30,
%     elseif D==50,load rosenbrock_M_D50,
%     else 
%         M=rot_matrix(D,c);
%         M=M.*(1+0.3.*normrnd(0,1,D,D));
%     end
    %----------------------------------------------------
end
x=x-repmat(o,ps,1);
% x=x*M; %---------------------------------------------------
fit=sum(x.^2,2);

% 	2.Shifted Schwefel's Problem 1.2
function f=schwefel_102(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load schwefel_102_data
   if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end

% 	3.Shifted Rotated High Conditioned Elliptic Function
function fit=high_cond_elliptic_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load high_cond_elliptic_rot_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    c=1;
    if D==2,load elliptic_M_D2,
    elseif D==10,load elliptic_M_D10,
    elseif D==30,load elliptic_M_D30,
    elseif D==50,load elliptic_M_D50,
%     elseif D==100,load elliptic_M_D100,
    else 
        A=normrnd(0,1,D,D);[M,r]=cGram_Schmidt(A);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M;
a=1e+6;
fit=0;
for i=1:D
fit=fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end

% 	4.Shifted Schwefel's Problem 1.2 with Noise in Fitness 
function f=schwefel_102_noise_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load schwefel_102_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
f=f.*(1+0.4.*abs(normrnd(0,1,ps,1)));

% 	5.Schwefel's Problem 2.6
function f=schwefel_206(x)%after Fletcher and Powell
global initial_flag
persistent A B o
[ps,D]=size(x);
if initial_flag==0
    initial_flag=1;
    load schwefel_206_data
    if length(o)>=D
         A=A(1:D,1:D);o=o(1:D);
    else
         o=-100+200*rand(1,D);
         A=round(-100+2*100.*rand(D,D));
         while det(A)==0
         A=round(-100+2*100.*rand(D,D));
         end
    end
    o(1:ceil(D/4))=-100;o(max(floor(0.75*D),1):D)=100;
    B=A*o';
end
for i=1:ps
f(i,1)=max(abs(A*(x(i,:)')-B));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Multimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	6.Shifted Rosenbrock's Function  ( 修改为 Rotated )
function f=rosenbrock_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load rosenbrock_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-90+180*rand(1,D);
    end
%     %------------------------------------------------
%     c=1;  %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
%     if D==2,load rosenbrock_M_D2,
%     elseif D==10,load rosenbrock_M_D10,
%     elseif D==30,load rosenbrock_M_D30,
%     elseif D==50,load rosenbrock_M_D50,
% %     elseif D==100,load rosenbrock_M_D100,
%     else 
%        M=rot_matrix(D,c);
%     end
%     %-------------------------------------------------    
%     o=o(1:D);
    initial_flag=1;
end
x=x-repmat(o,ps,1)+1;
% x=x*M; %-------------------------------------------------
f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);

% 	7.Shifted Rotated Griewank's Function
function f=griewank_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load griewank_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-600+0*rand(1,D);
    end
    c=3;   %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
    if D==2,load griewank_M_D2,  
    elseif D==10,load griewank_M_D10,  
    elseif D==30,load griewank_M_D30, 
    elseif D==50,load griewank_M_D50,  
%     elseif D==100,load griewank_M_D100,
    else 
        M=rot_matrix(D,c);
        M=M.*(1+0.3.*normrnd(0,1,D,D));
    end
   %----------------------------------------------------
    o=o(1:D);
    initial_flag=1;
end
x=x-repmat(o,ps,1); 
x=x*M;    %---------------------------------------------
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;

% 	8.Shifted Rotated Ackley's Function with Global Optimum on Bounds
function f=ackley_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load ackley_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-30+60*rand(1,D);
    end
    o(2.*[1:floor(D/2)]-1)=-32;
    c=100;%----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
    if D==2,load ackley_M_D2,
    elseif D==10,load ackley_M_D10,
    elseif D==30,load ackley_M_D30,
    elseif D==50,load ackley_M_D50,
%     elseif D==100,load ackley_M_D100,
    else 
       M=rot_matrix(D,c);
    end
    %-------------------------------------------------
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M; %-------------------------------------------------
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);

% 	9.Shifted Rastrign's Function
function f=rastrigin_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-5+10*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);

% 	10.Shifted Rotated Rastrign's Function 
function f=rastrigin_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-5+10*rand(1,D);
    end
    c=2;%----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
% %     ---------------------------------------------------
    if D==2,load rastrigin_M_D2,
    elseif D==10,load rastrigin_M_D10,
    elseif D==30,load rastrigin_M_D30,
    elseif D==50,load rastrigin_M_D50,
%     elseif D==100,load rastrigin_M_D100,
    else 
        M=rot_matrix(D,c);
    end
    %---------------------------------------------------
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M;%---------------------------------------------------
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);

% 	11.Shifted Rotated Weierstrass Function
function [f]=weierstrass_rot(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load weierstrass_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-0.5+0.5*rand(1,D);
    end
    c=5;%----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
    if D==2,load weierstrass_M_D2,
    elseif D==10,load weierstrass_M_D10,
    elseif D==30,load weierstrass_M_D30,
    elseif D==50,load weierstrass_M_D50,
%     elseif D==100,load weierstrass_M_D100,
    else 
        M=rot_matrix(D,c);
    end
    %---------------------------------------------------
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M; %---------------------------------------------------
x=x+0.5;
a = 0.5;%0<a<1
b = 3;
kmax = 20;
[ps,D]=size(x);

c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
c=-w(0.5,c1,c2);
f=0;
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+repmat(c*D,ps,1);

%--------------------------------

% 	12.Schwefel's Problem 2.13
function f=schwefel_213(x)%after Fletcher and Powell
global initial_flag
persistent a b A alpha
[ps,D]=size(x);
if initial_flag==0
    initial_flag=1;
    load schwefel_213_data
    if length(alpha)>=D
        alpha=alpha(1:D);a=a(1:D,1:D);b=b(1:D,1:D);
    else
        alpha=-3+6*rand(1,D);
        a=round(-100+200.*rand(D,D));
        b=round(-100+200.*rand(D,D));
    end
    alpha=repmat(alpha,D,1);
    A=sum(a.*sin(alpha)+b.*cos(alpha),2);
end

for i=1:ps
    xx=repmat(x(i,:),D,1);
    B=sum(a.*sin(xx)+b.*cos(xx),2);
    f(i,1)=sum((A-B).^2,1);
end

% 	13. Expanded Extended Griewank's plus Rosenbrock's Function (F8F2)
function fit=EF8F2_func(x)
%-3,1
global initial_flag
persistent  o 
[ps,D]=size(x);
if initial_flag==0
    load EF8F2_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-1+1*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1)+1;
fit=0;
for i=1:(D-1)
    fit=fit+F8F2(x(:,[i,i+1]));
end
    fit=fit+F8F2(x(:,[D,1]));
    
function f=F8F2(x)
f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
f=1+f2.^2./4000-cos(f2);

% ---------------------------------------------------------------  
% 	14. Expanded Rotated Extended Scaffer's F6 	
function f=E_ScafferF6_func(x)
global initial_flag
persistent  o M
fhd=str2func('ScafferF6');
[ps,D]=size(x);
if initial_flag==0
    load E_ScafferF6_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
    c=3;
    if D==2,load E_ScafferF6_M_D2,
    elseif D==10,load E_ScafferF6_M_D10,
    elseif D==30,load E_ScafferF6_M_D30,
    elseif D==50,load E_ScafferF6_M_D50,
%     elseif D==100,load E_ScafferF6_M_D100,
    else 
       M=rot_matrix(D,c);
    end
end
x=x-repmat(o,ps,1);
x=x*M;
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));

function f=ScafferF6(x)
f=0.5+(sin(sqrt(x(:,1).^2+x(:,2).^2)).^2-0.5)./(1+0.001*(x(:,1).^2+x(:,2).^2)).^2;
%---------------------------------------------------
%   15.Hybrid Composition Function 1
function fit=hybrid_func1(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func1_data % saved the predefined optima
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('frastrigin');
    func.f2=str2func('frastrigin');
    func.f3=str2func('fweierstrass');
    func.f4=str2func('fweierstrass');
    func.f5=str2func('fgriewank');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fackley');
    func.f8=str2func('fackley');
    func.f9=str2func('fsphere');
    func.f10=str2func('fsphere');
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num);
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    for i=1:fun_num
        eval(['M.M' int2str(i) '=diag(ones(1,D));']);
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%---------------------------------------------------------------------
%   16.Rotated Hybrid Composition Function 1	
function fit=hybrid_rot_func1(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func1_data % saved the predefined optima
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('frastrigin');
    func.f2=str2func('frastrigin');
    func.f3=str2func('fweierstrass');
    func.f4=str2func('fweierstrass');
    func.f5=str2func('fgriewank');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fackley');
    func.f8=str2func('fackley');
    func.f9=str2func('fsphere');
    func.f10=str2func('fsphere');
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num); 
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    c=[2,2,2,2,2,2,2,2,2,2,2];
    if D==2,load hybrid_func1_M_D2,
    elseif D==10,load hybrid_func1_M_D10,
    elseif D==30,load hybrid_func1_M_D30,
    elseif D==50,load hybrid_func1_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   17.	Rotated Hybrid Composition Function 1 with Noise in Fitness	
function fit=hybrid_rot_func1_noise(x)
[ps,D]=size(x);
fit=hybrid_rot_func1(x).*(1+0.2.*abs(normrnd(0,1,ps,1)));
%----------------------------------------------------------------
%   18.	Rotated Hybrid Composition Function 2
function fit=hybrid_rot_func2(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func2_data % saved the predefined optima
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    o(10,:)=0;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,load hybrid_func2_M_D2,
    elseif D==10,load hybrid_func2_M_D10,
    elseif D==30,load hybrid_func2_M_D30,
    elseif D==50,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   19.	Rotated Hybrid Composition Function 2 with a Narrow Basin for the Global Optimum
function fit=hybrid_rot_func2_narrow(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func2_data % saved the predefined optima
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    o(10,:)=0;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[0.1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[0.1*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,load hybrid_func2_M_D2,
    elseif D==10,load hybrid_func2_M_D10,
    elseif D==30,load hybrid_func2_M_D30,
    elseif D==50,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   20.	Rotated Hybrid Composition Function 2 with the Global Optimum on the Bounds	
function fit=hybrid_rot_func2_onbound(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func2_data % saved the predefined optima,
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    o(10,:)=0;
    o(1,2.*[1:floor(D/2)])=5;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,load hybrid_func2_M_D2,
    elseif D==10,load hybrid_func2_M_D10,
    elseif D==30,load hybrid_func2_M_D30,
    elseif D==50,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-------------------------------------------------
%    21.Rotated Hybrid Composition Function 3		
function fit=hybrid_rot_func3(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=ones(1,D);
    if D==2,load hybrid_func3_M_D2,
    elseif D==10,load hybrid_func3_M_D10,
    elseif D==30,load hybrid_func3_M_D30,
    elseif D==50,load hybrid_func3_M_D50,
    else 
        for i=1:fun_num
            A=normrnd(0,1,D,D);
            eval(['M.M' int2str(i) '=cGram_Schmidt(A));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-----------------------------------------
%   22.	Rotated Hybrid Composition Function 3 with High Condition Number Matrix
function fit=hybrid_rot_func3_highcond(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=[10 20 50 100 200 1000 2000 3000 4000 5000];
    if D==2,load hybrid_func3_HM_D2,
    elseif D==10,load hybrid_func3_HM_D10,
    elseif D==30,load hybrid_func3_HM_D30,
    elseif D==50,load hybrid_func3_HM_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-----------------------------------------
%   23.	Non-Continuous Rotated Hybrid Composition Function 3
function fit=hybrid_rot_func3_noncont(x)
global initial_flag
persistent  o 
[ps,D]=size(x);
if initial_flag==0
    load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
    o=o(1,1:D);
end
o=repmat(o,ps,1);
x=(abs(x-o)<0.5).*x+(abs(x-o)>=0.5).*(round(x.*2)./2);
fit=hybrid_rot_func3(x);
%-----------------------------------------
%   24.	Rotated Hybrid Composition Function 4	
function fit=hybrid_rot_func4(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    load hybrid_func4_data % saved the predefined optima, a 10*1000 matrix;
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fweierstrass');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('fEF8F2');
    func.f4=str2func('fackley');
    func.f5=str2func('frastrigin');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fE_ScafferF6_noncont');
    func.f8=str2func('frastrigin_noncont');
    func.f9=str2func('felliptic');
    func.f10=str2func('fsphere_noise');
    bias=((1:fun_num)-1).*100;
    sigma=[2,2,2,2,2,2,2,2,2,2];
    lamda=[10; 5/20; 1; 5/32; 1; 5/100 ; 5/50; 1; 5/100; 5/100; ];
    lamda=repmat(lamda,1,D);
    c=[100 50 30 10 5 5 4 3 2 2];
    if D==2,load hybrid_func4_M_D2,
    elseif D==10,load hybrid_func4_M_D10,
    elseif D==30,load hybrid_func4_M_D30,
    elseif D==50,load hybrid_func4_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------
function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M)
[ps,D]=size(x);
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2));
end

[tmp,tmpid]=sort(weight,2);
for i=1:ps
    weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
end
weight=weight./repmat(sum(weight,2),1,fun_num);

fit=0;
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1))*M.M' int2str(i) ');']);
    x1=5*ones(1,D);
    eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))*M.M' int2str(i) ');']);
    fit1=2000.*f./f1;
    fit=fit+weight(:,i).*(fit1+bias(i));
end
%-------------------------------------------------

%    46. frequency_modulated
function fit=frequency_modulated(x, o, A, M, a, alpha, b)
global initial_flag
theta=2*pi/100;
fit=0;
for t=0:100
   y_t=x(1)*sin(x(2)*t*theta+x(3)*sin(x(4)*t*theta+x(5)*sin(x(6)*t*theta)));
   y_0_t=1*sin(5*t*theta-1.5*sin(4.8*t*theta+2*sin(4.9*t*theta)));
   fit=fit+(y_t-y_0_t)^2;
end

%    47. shifted nonc_rastrigin_func
function fit=nonc_rastrigin_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-5+10*rand(1,D);
    end
    c=2;%----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
% %     ---------------------------------------------------
    if D==2,load rastrigin_M_D2,
    elseif D==10,load rastrigin_M_D10,
    elseif D==30,load rastrigin_M_D30,
    elseif D==50,load rastrigin_M_D50,
    else 
        M=rot_matrix(D,c);
    end
    %---------------------------------------------------
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
fit=sum(x.^2-10.*cos(2.*pi.*x)+10,2);


%    48. shifted and rotated  nonc_rastrigin_rot_func
function fit=nonc_rastrigin_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-5+10*rand(1,D);
    end
    c=2;%----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
% %     ---------------------------------------------------
    if D==2,load rastrigin_M_D2,
    elseif D==10,load rastrigin_M_D10,
    elseif D==30,load rastrigin_M_D30,
    elseif D==50,load rastrigin_M_D50,
    else 
        M=rot_matrix(D,c);
    end
    %---------------------------------------------------
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M;%---------------------------------------------------
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
fit=sum(x.^2-10.*cos(2.*pi.*x)+10,2);

%    49. shifted and rotated  salomon_func
function fit=Shifted_Rotated_salomon_func(x)
    global initial_flag
	persistent o M
    [ps,D]=size(x);
    if initial_flag==0
     load sphere_func_data   %  CEC2005中，只支持10D，30D,50D
%    load sphere_shift_func_data   %  CEC2008中，只支持100D，500D,1000D
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;    
    c=1;   %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
        if D==2,load elliptic_M_D2,  
        elseif D==10,load elliptic_M_D10,  
        elseif D==30,load elliptic_M_D30, 
        elseif D==50,load elliptic_M_D50,  
        elseif D==100,load elliptic_M_D100,
        else 
            M=rot_matrix(D,c);
            M=M.*(1+0.3.*normrnd(0,1,D,D));
        end
        %----------------------------------------------------
    end
    x=x-repmat(o,ps,1);
    x=x*M; %---------------------------------------------------
    f1=2*pi.*sqrt(sum(x.^2,2));
    f2=0.1.*sqrt(sum(x.^2,2))+1;
    fit=-cos(f1)+f2;%+fbias;
    
% 	50.Shifted Rotated Sphere Function 
function fit=Shifted_Rotated_sphere_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
%     load sphere_func_data   %  CEC2005中，只支持10D，30D,50D
    load sphere_shift_func_data   %  CEC2008中，只支持100D，500D,1000D
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
    
    c=1;   %----------- 若取消 Rotated ,则注释掉下列含“x=x*M;”的语句
    %---------------------------------------------------
    if D==2,load rosenbrock_M_D2,
    elseif D==10,load rosenbrock_M_D10,
    elseif D==30,load rosenbrock_M_D30,
    elseif D==50,load rosenbrock_M_D50,
    else 
        M=rot_matrix(D,c);
        M=M.*(1+0.3.*normrnd(0,1,D,D));
    end
    %----------------------------------------------------
end
x=x-repmat(o,ps,1);
x=x*M; %---------------------------------------------------
fit=sum(x.^2,2);

% 	51.sprd_spectrum_rad_pphase  % CEC2011-F7 
function fit=sprd_spectrum_rad_pphase(x)
    [ps,d]=size(x);
    hsum=[];
    var=2*d-1;
    for kk=1:2*var
      if rem(kk,2)~=0
         i=(kk+1)/2;
         hsum(kk)=0;
      for j=i:d    %fi(2i-1)X
         summ=0;
          for i1=(abs(2*i-j-1)+1):j
             summ=x(i1)+summ;
          end
         hsum(kk)= cos(summ)+ hsum(kk);
      end
      else 
        i=kk/2;
        hsum(kk)=0;
      for j=(i+1):d    %fi(2i)X
        summ=0;
         for i1=(abs(2*i-j)+1):j
           summ=x(i1)+summ;
         end
        hsum(kk)= cos(summ)+ hsum(kk);
      end
       hsum(kk)=hsum(kk)+0.5;
      end
    end
    
    fit=max(hsum) ;

% 	52.lennard_jones_potential_problem  % CEC2011-F2 
function fit=lennard_jones_potential_problem(x)
     [ps,d]=size(x);
      r=[];
      p=size(x);
      if rem(p(2),3)~=0
         disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
      end
       n=p(2)/3;
       x=reshape(x,3,n)';
       v=0;
       a=ones(n,n);
       b=repmat(2,n,n);
     for i=1:(n-1)
       for j=(i+1):n
        r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
        v=v+(a(i,j)/r(i,j)^12-b(i,j)/r(i,j)^6);
       end
     end
     fit=v;    
    
    
%basic functions

function f=fsphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D]=size(x);
f=sum(x.^2,2);
%--------------------------------
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
%--------------------------------
function f=fgriewank(x)
[ps,D]=size(x);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;
%--------------------------------
function f=fackley(x)
[ps,D]=size(x);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
%--------------------------------
function f=frastrigin(x)
[ps,D]=size(x);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function [f]=fweierstrass(x)
[ps,D]=size(x);
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c=-w(0.5,c1,c2);
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+c*D;

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum(c1 .* cos(c2.*x(:,k)));
end
%--------------------------------
function f=fE_ScafferF6(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%--------------------------------    
function f=fE_ScafferF6_noncont(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
    f=f+F8F2(x(:,[D,1]));

%--------------------------------
function f=fschwefel_102(x)
[ps,D]=size(x);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
%--------------------------------
function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
%--------------------------------

% classical Gram Schmid 
 function [q,r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 
% 
 [n,m] = size(A); 
 q = A;    
 for j=1:m
     for i=1:j-1 
         r(i,j) = q(:,j)'*q(:,i);
     end
     for i=1:j-1   
       q(:,j) = q(:,j) -  r(i,j)*q(:,i);
     end
     t =  norm(q(:,j),2 ) ;
     q(:,j) = q(:,j) / t ;
     r(j,j) = t  ;
 end 
 
function M=rot_matrix(D,c)
A=normrnd(0,1,D,D);
P=cGram_Schmidt(A);
A=normrnd(0,1,D,D);
Q=cGram_Schmidt(A);
u=rand(1,D);
D=c.^((u-min(u))./(max(u)-min(u)));
D=diag(D);
M=P*D*Q; 

function y = u(x,a,k,m)
y = zeros(length(x),1);
for kk = 1:length(x)
    if x(1,kk)>a
        y(kk)=k*((x(1,kk)-a)^m);
    elseif x(1,kk)<-a
        y(kk)=k*((-x(1,kk)-a)^m);
    else
        y(kk)=0;
    end	
end


% function y = w(x,c1,c2)
% y = zeros(length(x),1);
% for k = 1:length(x)
% 	y(k) = sum(c1 .* cos(c2.*x(:,k)));
% end



% if f<best_f
%     best_f=f;
% end
% best_keep=[best_keep,best_f];

