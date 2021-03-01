clear all;
close all;
clc;
T=2; %stage
A=2; %player 1's actions
B=2; %player 2's actions
k=3; %state of player 1
l=2; %state of player 2
lm=.9; %discounted value
load M.mat;  %payoff matrix
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
% p(1)=1; %initial probability p(k) when k=1
% p(2)=0; %initial probability p(k) when k=2
% p(3)=0; %initial probability p(k) when k=2
% q(1)=0; %initial probability q(l) when l=1
% q(2)=1; %initial probability q(l) when l=2


%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);


%Contraint R_{I_t}(a_t)=P_{a_{t-1},b_{t-1}}(k_{t-1},k_t)R_{I_{t-1}}(a_t-1)
Aeq=zeros(sum(n_is1),sum(n_is1)*A);
beq=zeros(sum(n_is1),1);

row_index=0;
for t=1:T
    for i=1:length(is1{t})
        row_index=row_index+1;
        for a=1:A  %As there is summation a_t in R_{I_{t}}
            I=is1{t}(i,:); %for selecting each row from the information set of player 1
            [col_index_RIt] = RIt_col_index_P1(t,I,A,B,k,a,n_is1); %finding the column index for that particular information set row
            Aeq(row_index,col_index_RIt)=1;
        end
    kt=I(end);
        if t==1
           beq(row_index,1)=p(kt); %As at t=1, P_{a_0,b_0}(k_0,k_1)=1 and R_{I_0}(a_0)=Pr(k1)
        else
            Ipre=is1{t}(i,1:(length(is1{t}(i,:))-3));
            aprev=is1{t}(i,(length(is1{t}(i,:))-2));
            bprev=is1{t}(i,(length(is1{t}(i,:))-1));
            [col_index_RItprev] = RIt_col_index_P1(t-1,Ipre,A,B,k,aprev,n_is1);
            ktpre=I(end-3);
            Aeq(row_index,col_index_RItprev)=-P{aprev,bprev}(ktpre,kt);
        end     
    end
end  


%for inequality contraints
[kset,n_kset] = Kset(T,k);
[J,n_J]=info_J(T,A,B,l);

Ain=zeros(sum(n_is2)*B,(sum(n_is1)*A+sum(n_is2)));
bin=zeros(sum(n_is2)*B,1);

row_index=0;
for t=1:T
for j=1: length (J{t}) %for different J_t and different b_t there will be different row
    for b=1:B
        row_index=row_index+1;
        for kn=1:length(kset{t})
            for a=1:A
                [I] = construct_Is(j,kn,t,kset,J);
                [col_index] = RIt_col_index_P1(t,I,A,B,k,a,n_is1);
                Ain(row_index,col_index)=lm^(t-1)*G{I(3*t-2),(J{t}(j,(3*t-2)))}(a,b); %R_{I_s} coefficient
            end
        end
        %To find the term for U_{J_{s+1}}
        if t<T  %As for U_{J_{N+1}}=0
            for a1=1:A
                for lplus=1:l  
                    Jplus=[J{t}(j,:) a1 b lplus];
                    [col_index_Jplus] = J_col_index(n_is1,A,B,l,Jplus,n_is2,t+1);
                     Ain(row_index,col_index_Jplus)=Q{a1,b}(J{t}(j,end),lplus);%The coefficient for U_{J_{s+1}}
                end
            end      
        end
        %To find the term for U_{J_{s}}
        Jpre=J{t}(j,:);
        [col_index_Jpre] = J_col_index(n_is1,A,B,l,Jpre,n_is2,t);
         Ain(row_index,col_index_Jpre)=-1;
    end
end
end 

%Construck the objective function
ql=zeros(1,l);
n=1;
for l_initial=1:l
    ql(1,n)=q(l_initial);
    n=n+1;
end
f=[zeros(1,sum(n_is1)*A) ql zeros(1,sum(n_is2(2:end)))]; 

%Rearrange every coefficient to use in linprog
Ain1=-Ain; %As linprog use <=
bin1=-bin; %As linprog use <=
f1=-[zeros(1,size(Ain1,2)-size(f,2)) f ]; %as linprog works for only minnimize obj function thats why -ve
Aeq1=padarray(Aeq,[0, (size(Ain1,2)-size(Aeq,2))],0,'post');
lb=[zeros(sum(n_is1)*A,1);-Inf((size(Ain1,1)-sum(n_is1)*A),1)];
ub=+Inf;

[x,v]=linprog(f1,Ain1,bin1,Aeq1,beq,lb,ub)

%Finding the optimal strategy from the realization plan we got from previous linprog result
%sigma=RIT/(P_{a_{t-1},b_{t-1}}(k_{t-1},k_t)*RI_{t-1})

sigma=zeros(1,sum(n_is1));
sigma_col=0;
for t=1:T
    for i=1:length(is1{t})
        I=is1{t}(i,:); %for selecting each row from the information set of player 1
        sigma_col=sigma_col+1;
        for a=1:A  %As there is summation a_t in R_{I_{t}}
            [row_index_RIt] = RIt_col_index_P1(t,I,A,B,k,a,n_is1); %finding the column index for that particular information set row. The number of column index in RIt is the number of row index in optimal realization plan x
            RIt=x(row_index_RIt,1);        
            kt=I(end);
            if t==1
               RItprev=p(kt);
               pvalue=1;
            else
                Ipre=is1{t}(i,1:(length(is1{t}(i,:))-3));
                aprev=is1{t}(i,(length(is1{t}(i,:))-2));
                bprev=is1{t}(i,(length(is1{t}(i,:))-1));
                [row_index_RItprev] = RIt_col_index_P1(t-1,Ipre,A,B,k,aprev,n_is1);
                ktpre=I(end-3);
                RItprev=x(row_index_RItprev,1);
                pvalue=P{aprev,bprev}(ktpre,kt);
            end
            sigma(a,sigma_col)=RIt/(RItprev*pvalue);
        end
    end
end
save sigma.mat sigma;

%to find ontimal nu
%to find UJ1
for l_present=1:l
    t=1;
    info=l_present;
    [col_index_UJ1] = J_col_index(n_is1,A,B,l,info,n_is2,t);
    UJ1=x(col_index_UJ1,1);
    nu(l_present)=[UJ1];
end

save nu.mat 'nu'