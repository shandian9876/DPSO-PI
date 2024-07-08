
function [exemplar, prob_1, prob_2]=Generate(iter,me,Queue_1,prob_1,pos_1, Queue_2,prob_2,pos_2,   Xmax,Xmin)
    N=size(Queue_1,1);
    D=size(pos_1,2);   
    pc=0.9-0.4*sin((iter/me)*pi); 
    pm=0.02; 
    p_1=prob_1;
    for i=2:N
       p_1(i)=p_1(i)+p_1(i-1); 
    end
    p_2=prob_2;
    for i=2:N
       p_2(i)=p_2(i)+p_2(i-1);
    end
    for j=1:D  
        for i=1:N
           if rand<=p_1(i)
               index=i;
               x_1=pos_1(index,:);
               break;
           end
        end    
        for i=1:N
           if rand<=p_2(i)
               index=i;
               x_2=pos_2(index,:);
               break;
           end
        end
        rd=rand;
        rd=rd<pc;   
        exemplar(1,j)=rd.*x_1(1,j)+(1-rd).*x_2(1,j);  
    end
    if rand >= iter/me
    p_3=prob_2;
    for i=2:N
       p_3(i)=p_3(i)+p_3(i-1); 
    end
        for j=1:D
          for i=1:N
           if rand<=p_3(i)
               index=i; 
               x_3=pos_1(index,:);
               break;
           end
          end
        rd=rand;
        rd=rd<pc;  
        exemplar(1,j)=rd.*exemplar(1,j)+(1-rd).*x_3(1,j);  
        end
    end
    for j=1:D
       if rand<pm
           exemplar(1,j)=Xmin(1,j)+(Xmax(1,j)-Xmin(1,j)).*rand;  
       end
    end
end