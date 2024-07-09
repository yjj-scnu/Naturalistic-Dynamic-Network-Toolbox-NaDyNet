function average_dwell_time=sf_ave_dwell_time(labels,K,TR)
T=length(labels);
for k=1:K   
    n_cluster=0;
    
       
    for j=2:T
        if labels(j)==k&&labels(j-1)~=k         %the adjacent labels is same???
            n_cluster=n_cluster+1;
        end
        
    end
   
        
   states_location=find(labels==k);
   if n_cluster
       average_dwell_time(k)=(length(states_location)*TR)/n_cluster;
   else
       average_dwell_time(k)=0;
   end
end


