function exemplar=DE_based_exemplar(N,D,Pbest,Xmax,Xmin)
    f=0.7;
    OP = KNN(Pbest);
    for i = 1:N
         neighbour_index=get_neighbor(i,N,OP);
        
%        Mut_Vec1=Pbest(i,:)+ f.*(Pbest(OP(ind-1),:) - Pbest(OP(ind+1),:));%变异
       Mut_Vec1=Pbest(i,:)+f.*(Pbest(neighbour_index(1),:)-Pbest(neighbour_index(3),:));
        pc=0.1;
        pro=rand(1,D)>pc;
        Mut_vec=pro.*Mut_Vec1+(1-pro).*Pbest(i,:); %交叉
        Mut_vec(Mut_vec>Xmax) = Xmax;
        Mut_vec(Mut_vec<Xmin) = Xmin;
        exemplar(i,:)=Mut_vec;
    end
end