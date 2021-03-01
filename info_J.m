function [is2 n_is2]=info_J(T,A,B,l)
index=0;
for lj=1:l
    index=index+1;
    is2{1}(index,:)=[lj];
end
n_is2(1)=index;
for t=2:T
    n_row=size(is2{t-1},1);
         index=0;
         for row=1:n_row % for every history information set
            for a=1:A; % new a
               for b=1:B; % new b
                   for lj=1:l %new info about state
                index=index+1;
                is2{t}(index,:)=[is2{t-1}(row,:) a b lj];
                   end
               end
            end
          end
    n_is2(t)=index; %number of row in each stage
    is2{t};  
end