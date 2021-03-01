function [b] = choose_action_randomly_p2(l_present, B)
    %setting a probability of choosing action according to l_present
    tau=[0.9 0.75; 0.1 0.25];
    %Choose player 1's action a_1 accoriding to the strategy sigma
    rd=rand(1);
    b=1;
    while sum(tau(1:b,l_present))<rd
        b=b+1;
        if b>B %For some reason we are getting a>A here. To solve this problem we have used a-1 and break
           b=b-1;
           break
        end
    end
end


