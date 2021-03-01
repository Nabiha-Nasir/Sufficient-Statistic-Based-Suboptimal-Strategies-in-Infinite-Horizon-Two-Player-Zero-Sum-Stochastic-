clear all;
close all;
clc;
T=2; %stage
A=2; %player 1's actions
B=2; %player 2's actions
k=3; %state of player 1
l=2; %state of player 2
lm=.9; %discounted value
load M.mat;     %payoff matrix
G=M;
P{1,1}=[.8 .1 .1; .1 .4 .5; .2 .7 .1 ]; % Payoff matrix if a=1 and b=1
P{1,2}=[.4 .5 .1; .2 .3 .5; .4 .4 .2 ]; % Payoff matrix if a=1 and b=2
P{2,1}=[.2 .2 .6; .5 .2 .3; .2 .2 .6]; % Payoff matrix if a=2 and b=1
P{2,2}=[.3 .3 .4; .1 .8 .1; .1 .1 .8 ]; % Payoff matrix if a=2 and b=2
Q{1,1}=[.8 .2;.5 .5]; % Payoff matrix if a=1 and b=1
Q{1,2}=[.2 .8;.1 .9]; % Payoff matrix if a=1 and b=2
Q{2,1}=[.6 .4;.5 .5]; % Payoff matrix if a=2 and b=1
Q{2,2}=[.7 .3; .1 .9]; % Payoff matrix if a=2 and b=2
save P.mat P
save Q.mat Q
% P{1,1}=[.8 .1 .1; .1 .4 .5; .2 .7 .1 ]; % Payoff matrix if a=1 and b=1
% P{1,2}=[.5 .5 0; 0 .5 .5; .4 .4 .2 ]; % Payoff matrix if a=1 and b=2
% P{2,1}=[.2 .2 .6; .5 .5 0; 0 0 1 ]; % Payoff matrix if a=2 and b=1
% P{2,2}=[.3 .3 .4; 0 1 0; .1 .1 .8 ]; % Payoff matrix if a=2 and b=2
% Q{1,1}=[1 0;.5 .5]; % Payoff matrix if a=1 and b=1
% Q{1,2}=[.2 .8;.1 .9]; % Payoff matrix if a=1 and b=2
% Q{2,1}=[.6 .4;.5 .5]; % Payoff matrix if a=2 and b=1
% Q{2,2}=[.7 .3; 0 1]; % Payoff matrix if a=2 and b=2
p(1)=0.5; %initial probability p(k) when k=1
p(2)=0.3; %initial probability p(k) when k=2
p(3)=0.2; %initial probability p(k) when k=2
q(1)=0.5; %initial probability q(l) when l=1
q(2)=0.5; %initial probability q(l) when l=2


%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);


%Contraint S_{J_t}(b_t)=Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)S_{J_{t-1}}(b_t-1)
Aeq=zeros(sum(n_is2),sum(n_is2)*B);
beq=zeros(sum(n_is2),1);


row_index=0;
for t=1:T
    for j=1:length(is2{t})
        row_index=row_index+1;
        for b=1:B  %As there is summation b_t in S_{J_{t}}
            J=is2{t}(j,:); %for selecting each row from the information set of player 2
            [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
            Aeq(row_index,col_index_SJt)=1;
        end
        lt=J(end);
        if t==1
           beq(row_index,1)=q(lt); %As at t=1, Q_{a_0,b_0}(l_0,l_1)=1 and S_{J_0}(b_0)=Pr(l1)
        else
            Jpre=is2{t}(j,1:(length(is2{t}(j,:))-3));
            aprev=is2{t}(j,(length(is2{t}(j,:))-2));
            bprev=is2{t}(j,(length(is2{t}(j,:))-1));
            [col_index_SJtprev] = SJt_col_index_P2(t-1,Jpre,A,B,l,bprev,n_is2);
            ltpre=J(end-3);
            Aeq(row_index,col_index_SJtprev)=-Q{aprev,bprev}(ltpre,lt);
        end
    end
end

    
%for inequality contraints
[lset,n_lset] = Lset(T,l);
[I,n_I]=info_I(T,A,B,k);

Ain=zeros(sum(n_is1)*A,(sum(n_is2)*B+sum(n_is1)));
bin=zeros(sum(n_is1)*A,1);

row_index=0;
for t=1:T
for i=1: length (I{t})
    for a=1:A
        row_index=row_index+1;
        for ln=1:length(lset{t})
            for b=1:B
                [J] = construct_Is(i,ln,t,lset,I); %As I and J are similar so we can use same function to construc both
                [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
                Ain(row_index,col_index_SJt)=lm^(t-1)*G{(I{t}(i,(3*t-2))),J(3*t-2)}(a,b);
            end
        end
        if t<T
            for b1=1:B
                for kplus=1:k  
                    Iplus=[I{t}(i,:) a b1 kplus]; %problem in this line. do I need to use different row for different b1?
                    [col_index_Isplus] = I_col_index(n_is1,A,B,k,Iplus,n_is2,t+1);
                     Ain(row_index,col_index_Isplus)=P{a,b1}(I{t}(i,end),kplus);
                end
            end     
        end
        Ipre=I{t}(i,:);
        [col_index_Ipre] = I_col_index(n_is1,A,B,k,Ipre,n_is2,t);
         Ain(row_index,col_index_Ipre)=-1;
    end
end
end

%Construct the objective function
pk=zeros(1,l);
n=1;
for k_initial=1:k
    pk(1,n)=p(k_initial);
    n=n+1;
end
f=[zeros(1,sum(n_is2)*B) pk zeros(1,sum(n_is1(2:end)))]; 

%Rearrange every coefficient to use in linprog
f1=[zeros(1,size(Ain,2)-size(f,2)) f ]; 
Aeq1=padarray(Aeq,[0, (size(Ain,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
lb=[zeros(sum(n_is2)*B,1);-Inf((size(Ain,1)-sum(n_is2)*B),1)];
ub=+Inf;

[x,v]=linprog(f1,Ain,bin,Aeq1,beq,lb,ub)

%Finding the optimal strategy from the realization plan we got from previous linprog result 
%tau=SJT/(Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)*SJ_{t-1})

tau=zeros(1,sum(n_is2));
tau_col=0;

row_index=0;
for t=1:T
    for j=1:length(is2{t})
        tau_col=tau_col+1;
        J=is2{t}(j,:); %for selecting each row from the information set of player 2
        for b=1:B
            [row_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
            SJt=x(row_index_SJt,1);
        lt=J(end);
        if t==1
            SJtprev=q(lt);
               qvalue=1;
        else
            Jpre=is2{t}(j,1:(length(is2{t}(j,:))-3));
            aprev=is2{t}(j,(length(is2{t}(j,:))-2));
            bprev=is2{t}(j,(length(is2{t}(j,:))-1));
            [row_index_SJtprev] = SJt_col_index_P2(t-1,Jpre,A,B,l,bprev,n_is2);
            ltpre=J(end-3);
            SJtprev=x(row_index_SJtprev,1);
            qvalue=Q{aprev,bprev}(ltpre,lt);
        end
          tau(b,tau_col)=SJt/(SJtprev*qvalue);
       end
    end
end
save tau.mat tau;

%to find mu
%to find ZI1
for k_present=1:k
    t=1;
    info=k_present;
    [col_index_ZI1] = I_col_index(n_is1,A,B,k,info,n_is2,t);
    ZI1=x(col_index_ZI1,1);
    mu(k_present)=[ZI1];
end

save mu.mat 'mu'