function [I] = construct_Is(j,kn,t,kset,J)
if t==1
    I= kset{t}(kn,:);
else
    n=1;
    Ip{t}(j,:)= J{t}(j,:); %As I need to replace l1, l2...ln in J_t to construct I. Thats why giving it a different name so that J_t can be unchanged
      for m=1:3:length(J{t}(j,:))
          Ip{t}(j,m)= kset{t}(kn,n); %replaceing l by k
          n=n+1;
          I=Ip{t}(j,:);
      end   
end
end

