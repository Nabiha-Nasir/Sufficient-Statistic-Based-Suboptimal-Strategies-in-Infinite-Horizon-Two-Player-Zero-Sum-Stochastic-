function [is1 n_is1]=info_I(T,A,B,k)
index=0;
for ki=1:k
    index=index+1;
    is1{1}(index,:)=[ki];
end
n_is1(1)=index;
for t=2:T
    n_row=size(is1{t-1},1);
         index=0;
         for row=1:n_row % for every history information set
            for a=1:A % new a
               for b=1:B % new b
                   for ki=1:k %new info about state
                index=index+1;
                is1{t}(index,:)=[is1{t-1}(row,:) a b ki];
                   end
               end
            end
          end
    n_is1(t)=index; %number of row in each stage
    is1{t};  
end

