function index=get_neighbor(i,N,indFitness)

    ind = find(i==indFitness);
    if ind==1 
        index=[indFitness(N),indFitness(ind),indFitness(ind+1)];
    elseif ind==N
        index=[indFitness(ind-1),indFitness(ind),indFitness(1)];
    else
        index=[indFitness(ind-1),indFitness(ind),indFitness(ind+1)];
    end
end

% function index=get_neighbor(i,N)
% if i==1
%     index=[N,1,2];
% elseif i==N
%     index=[i-1,N,1];
% else
%     index=[i-1,i,i+1];
% end