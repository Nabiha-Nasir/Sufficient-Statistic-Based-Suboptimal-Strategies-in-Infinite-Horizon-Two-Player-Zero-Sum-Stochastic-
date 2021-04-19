clear all;
close all;
clc;
warning('off')
%%% both players are using receding horizon to compute the optimal strategy part by part
%% initialization
load M.mat;%payoff matrix
G=M;
load P.mat %Player 1's transition matrix
load Q.mat %Player 2's transition matrix
N=36; %Game horizon
T=2; %window size. So in window by window it will be 1:2, 3:4
A=2; %Number of player 1's actions
B=2; %Number of player 2's actions
k=3; %Number of states of player 1
l=2; %Number of states of player 2
lm=0.3; %discounted value
p=[0.5 0.3 0.2]; %player 1's initial probability for state
q=[0.5 0.5];     %player 2's initial probability for state

total_loop=500; %number of experiemnts to get the average payoff
for_avg=0;
payoff_detail_all=[];

for loop_count=1:total_loop 
   
    % getting the first window strategy and initial vector payoff from primal game
    [sigma,nu] = fn_primal_game_p1(T,A,B,k,l,lm,P,Q,p,q,G);
    [tau,mu] = fn_primal_game_p2(T,A,B,k,l,lm,P,Q,p,q,G);
    p_present=p;
    q_present=q;
    
    payoff_detail=[];
    
    %% For N=1
    X_star=sigma(:,1:k); % This is the optimal strategy for player 1 in first stage
    Y_star=tau(:,1:l); %This is the optimal strategy for player 2 in first stage
    
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
    [col_index_tau] = tau_col_index_new(A,B,l,l_present,Hb,1); %We set t as 1. For t=1 and specific Hb this is my Y*
    [b] = choose_action(tau,B,col_index_tau);
    action_p2=[action_p2 b];

    %Payoff for stage 1
    payoff=G{k_present,l_present}(a,b)*lm^(1-1); %here we set the t in lm^(1-1) as 1. As it is the first stage
    payoff_detail=[payoff_detail payoff];
    
    %% for N=2:N
    stage_counter=2; %here n is just a operator to count the loop 
    while stage_counter<=N
        
        % to find mu
        Y_star(isnan(Y_star))=0; % to convert the NAN values to zero
        Y_star=abs(Y_star);
        %to update mu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [beta_vector] = fn_dual_game_2nd_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu,Y_star);
        %to choose mu corresponding to the previous actions
        col_index_beta=(a-1)*B*k+(b-1)*k;
        mu=beta_vector(1,col_index_beta+1:col_index_beta+k);

        %to find q+
        [q_present] = update_belief_li_qplus(l,q_present,Y_star,a,b,Q);

        %to find nu
        X_star(isnan(X_star))=0; % to convert the NAN values to zero
        X_star=abs(X_star);
        %to update nu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [alpha_vector] = fn_dual_game_2nd_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu,X_star);
        %to choose nu corresponding to the previous actions
        col_index_alpha=(a-1)*B*l+(b-1)*l;
        nu=alpha_vector(1,col_index_alpha+1:col_index_alpha+l);

        %to find p+
        [p_present] = update_belief_li_pplus(k,p_present, X_star,a,b,P);

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

        %Choose player 2's action accoriding to the strategy updated tau
        [col_index_tau_plus] = tau_col_index_new(A,B,l,l_present,Hb,1); %For t=1 and specific Hb this is my Y*. In this case no matter what the Hb is, it is not using it as t=1.
        [b] = choose_action(tau_plus,B,col_index_tau_plus);
        action_p2=[action_p2 b]; 

        %to find new \sigma*
        sigma_plus =fn_full_dual_game_value_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu);
        X_star=sigma_plus(:,1:k); %this means all rows and col 1: to k. here t=1

        %Choose player 1's action accoriding to the strategy updated sigma
        [col_index_sigma_plus] = sigma_col_index_new(A,B,k,k_present,Ha,1); %For t=1 and specific Ha this is my X*. In this case no matter what the Ha is, it is not using it as t=1.
        [a] = choose_action(sigma_plus,A,col_index_sigma_plus);
        action_p1=[action_p1 a]; 

        %Compute the payoff of every stage
        payoff=payoff+G{k_present,l_present}(a,b)*lm^(stage_counter-1);
        payoff_detail=[payoff_detail G{k_present,l_present}(a,b)*lm^(stage_counter-1)];
        stage_counter=stage_counter+1;
    end
    payoff_detail_all=[payoff_detail_all; payoff_detail];
    for_avg=for_avg+payoff;
end

payoff;
average_payoff=for_avg/total_loop
