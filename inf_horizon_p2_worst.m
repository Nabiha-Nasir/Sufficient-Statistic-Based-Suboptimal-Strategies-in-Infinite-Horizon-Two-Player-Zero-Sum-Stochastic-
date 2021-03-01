clear all;
close all;
clc;
%initialization
N=36; %stage. 
warning('off')
%we will roll for n=2. So the rolling horizon will be 1:2, 3:4,....19:20
T=2; %stage
A=2; %player 1's actions
B=2; %player 2's actions
k=3; %state of player 1
l=2; %state of player 2
lm=.9; %discounted value
load M.mat; %payoff matrix
G=M;
load P.mat;
load Q.mat;
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

%Before doing this test make sure the files from where we are loading mu,nu,sigma,tau are also for same T

%setting the number of runs
total_loop=500;

for_avg=0;
payoff_detail_all=[];


for M=1:total_loop %this loop is to run same experiment many times to get the avg value
    
    %loading data from primal_game_value_P1
    load mu;
    mu=-mu;
    load sigma.mat; %sigma is the strategy of player 1. Its each row is for each action (a) and each column is for each information set

    %loading data from primal_game_value_P2
    load nu;
    nu=-nu;
    load tau.mat;   %tau is the strategy of player 2. Its each row is for each action (b) and each column is for each information set
 
    p_present=p;
    q_present=q;
    
    payoff_detail=[];
    
    %% For N=1
    X_star=sigma(:,1:k); %this means all rows and col 1: to k. This is the optimal strategy for player 1 in t=1
    Y_star=tau(:,1:l); %this means all rows and col 1: to L. This is the optimal strategy for player 2 in t=1
    
    %Choose 1st stage's state of player 1: k_1
    [k_present] = choose_state(p,k);

    %Choose 1st stage's state of player 2: l_1
    [l_present] = choose_state(q,l);

    %updating history information
    Ha=[k_present];
    Hb=[l_present]; 
    action_p1=[];
    action_p2=[];

    %Choose player 1's action a_1 accoriding to the strategy sigma
    [col_index_sigma] = sigma_col_index_new(A,B,k,k_present,Ha,1); %We set t as 1. For t=1 and specific Ha this is my X*
    [a] = choose_action(sigma,A,col_index_sigma);
    action_p1=[action_p1 a];

    %Choose player 2's action b_1 accoriding to the strategy tau
    %[col_index_tau] = tau_col_index_new(A,B,l,l_present,Hb,1); %We set t as 1. For t=1 and specific Hb this is my Y*
    [b] = choose_action_randomly_p2(l_present, B);
    action_p2=[action_p2 b];

    %Payoff for stage 1
    payoff=G{k_present,l_present}(a,b)*lm^(1-1); %here we set the t in lm^(1-1) as 1. As it is the first stage
    payoff_detail=[payoff_detail payoff];

    
    %% for N=2:N
    stage_counter=2; %here n is just a operator to count the loop 
    while stage_counter<=N
        Y_star(isnan(Y_star))=0; % to convert the NAN values to zero
        Y_star=abs(Y_star);
    %to update mu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [beta_vector] = fn_dual_game_2nd_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu,Y_star);
    %to choose mu corresponding to the previous actions
        col_index_beta=(a-1)*B*k+(b-1)*k;
        mu=beta_vector(1,col_index_beta+1:col_index_beta+k);

    %to find q+
        [q_present] = update_belief_li_qplus(l,q_present,Y_star,a,b,Q);
%         if sum(q_present)<0.8
%             q_present
%         end

        X_star(isnan(X_star))=0; % to convert the NAN values to zero
        X_star=abs(X_star);
        if sum(X_star,'all')<.9
            stage_counter
            X_star
        end
    %to update nu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [alpha_vector] = fn_dual_game_2nd_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu,X_star);
    %to choose nu corresponding to the previous actions
        col_index_alpha=(a-1)*B*l+(b-1)*l;
        nu=alpha_vector(1,col_index_alpha+1:col_index_alpha+l);

    %to find p+
        [p_present] = update_belief_li_pplus(k,p_present, X_star,a,b,P)
        if sum(p_present)<0.8
            stage_counter
            p_present
        end

    % update state
    %for k
        k_prob=P{a,b}(k_present,:);
        [k_present] = choose_state(k_prob,k);

    %for l
        l_prob=Q{a,b}(l_present,:);
        [l_present] = choose_state(l_prob,l);

    % update history
          Ha=[Ha a b k_present];
          Hb=[Hb a b l_present];

    %to find new \tau*
        tau_plus = fn_full_dual_game_value_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu);
        Y_star=tau_plus(:,1:l); %this means all rows and col 1: to L. here t=1
        if sum(Y_star,'all')<1.5
            Y_star
        end

    %Choose player 2's action accoriding to the strategy updated tau
        %[col_index_tau_plus] = tau_col_index_new(A,B,l,l_present,Hb,1); %For t=1 and specific Hb this is my Y*. In this case no matter what the Hb is, it is not using it as t=1.
        [b] = choose_action_randomly_p2(l_present, B);
        action_p2=[action_p2 b]; 


    %to find new \sigma*
        sigma_plus =fn_full_dual_game_value_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu);
        X_star=sigma_plus(:,1:k); %this means all rows and col 1: to k. here t=1
        if sum(X_star,'all')<2
            X_star
        end

    %Choose player 1's action accoriding to the strategy updated sigma
        [col_index_sigma_plus] = sigma_col_index_new(A,B,k,k_present,Ha,1); %For t=1 and specific Ha this is my X*. In this case no matter what the Ha is, it is not using it as t=1.
        [a] = choose_action(sigma_plus,A,col_index_sigma_plus);
        action_p1=[action_p1 a]; 

    %Compute the payoff of every stage
        payoff=payoff+G{k_present,l_present}(a,b)*lm^(stage_counter-1);
        payoff_detail=[payoff_detail G{k_present,l_present}(a,b)*lm^(stage_counter-1)];
        stage_counter=stage_counter+1;
    end
for_avg=for_avg+payoff;
payoff_detail_all=[payoff_detail_all; payoff_detail];
end

payoff;
average_payoff=for_avg/total_loop
