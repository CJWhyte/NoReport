%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Foundational no-report simulation -------------%

% This script qualitatively reproduces the simulation results reported in 
% Whyte, Hohwy, and Smith (2021) using a model that also generates 
% "verbal reports" of the kind reported in Whyte and Smith (2020). 

% Code written by Christopher Whyte and Ryan Smith.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
rng('shuffle')

dbstop if error


%% task settings

% set to 1 for report condition, and 2 for no-report condition
report = 1; 

for cons = 1:2

    if cons == 1
        evi = 0.5; %conscious condition
    elseif cons == 2
        evi = 0.05; %unconscious condition
    end 
    
    %% Level 1
    %==========================================================================

    % prior beliefs about initial states
    %--------------------------------------------------------------------------

    D{1} = [1 1 1]';% Gabor stim {left, right, blank}
    D{2} = [1 0 0]';% saleince maps {null, ~maintain, maintain}
    D{3} = [1 0 0 0 0 0 0]';% language processing {silent, "I", "saw" "a", "Gabor", "didn't", "anything"}

    % probabilistic mapping from hidden states to outcomes: A
    %--------------------------------------------------------------------------

    %outcome modality 1: stimulus

    for i = 1:3
        for  j = 1:7
            A{1}(:,:,i,j) = eye(3,3);
        end 
    end 


    
    maintance = .9; %lower first level precision to allow WM maintance
    %seperate generative process from generative model
    a = A;

    %General decrease of stimulus strength - null
    a{1}(:,:,1) = spm_softmax(evi*log(a{1}(:,:,1)+exp(-4)));

    %WM attention modulation: reduction of first level precision
    a{1}(:,:,3) = spm_softmax(maintance*log(a{1}(:,:,3)+exp(-4)));


    %multiplying by 64 prevents learning by making concentration parameteres very high).
    a{1}= a{1}*64;

    % Transitions between states: B
    %--------------------------------------------------------------------------

    B{1}= eye(3,3);
    B{2} = eye(3,3);
    B{3} = eye(7,7);

    % MDP Structure
    %--------------------------------------------------------------------------
    mdp.T = 1;                      % number of updates
    mdp.A = A;                      % observation model
    mdp.B = B;                      % transition probabilities
    mdp.D = D;                      % prior over initial states
    mdp.a = a;
    mdp.erp = 1;

    mdp.Aname = {'Stimulus'};
    mdp.Bname = {'Stimulus', 'Attention', 'Langauge'};

    clear a A B D

    MDP = spm_MDP_check(mdp);

    clear mdp


    %% Level 2 (slower semantic timescale)
    %==========================================================================

    % prior beliefs about initial states (in terms of counts_: D and d
    %--------------------------------------------------------------------------
    D{1} = [1 1 1]'; % Sequence type: {left, right, blank}
    D{2} = [1 0 0 0 0 0 0 0 0]'; % time in trial: {blank, stimulus_1, stimulus_2, delay, delay, report}
    D{3} = [1 0 0]'; % Goal: {null, ~maintain, maintain}
    D{4} = [1 0 0]'; % Report: {null, unseen, seen} or {null, left, right} 
    d = D;

    % probabilistic mapping from hidden states to outcomes: A
    %--------------------------------------------------------------------------

    % outcomes: A{1} stim (3), A{2} metabolic demand (2), A{2} Report
    % Feedback (3) A{5} = report(2)

    %--- Stimulus
    for i = 1:9
        for j = 1:3
            for k = 1:3
                A{1}(:,:,i,j,k) = [0 0 0;%left
                                   0 0 0;%right
                                   1 1 1];%blank
            end 
        end
    end


    %present stimulus on trial 2 & 3
    for i = 2:3
          for j = 1:3
              for k = 1:3
                    A{1}(:,:,i,j,k) = [1 0 0;%left
                                       0 1 0;%right
                                       0 0 1];%blank
              end 
          end
    end


    %maintain stimulus after presentation
    for i = 4:8
          for j = 3
              for k = 1:3
                    A{1}(:,:,i,j,k) = [1 0 0;%left
                                       0 1 0;%right
                                       0 0 1];%blank
              end 
          end
    end


    %--- Goal
    for i = 1:9
          for j = 1:3
              for k = 1:3
                    A{2}(:,:,i,j,k) = [1 1 1;% null
                                       0 0 0;% low metabolic demand
                                       0 0 0];% high metabolic demand
              end 
          end
    end

    % ~ maintain
    for i = 4:8
          for j = 2
              for k = 1:3
                    A{2}(:,:,i,j,k) = [0 0 0;% null
                                       1 1 1;% low metabolic demand
                                       0 0 0];% high metabolic demand
              end 
          end
    end


    % maintain
    for i = 4:8
          for j = 3 
              for k = 1:3
                    A{2}(:,:,i,j,k) = [0 0 0;% null
                                       0 0 0;% low meabolic demand
                                       1 1 1];% high metabolic demand
              end 
          end
    end



    %--- Report
    for i = 1:9
        for j = 1:3 %"no goal"
            for k = 1:3
                  A{3}(:,:,i,j,k) = [1 1 1;%null
                                     0 0 0;%correct
                                     0 0 0];%incorrect
            end 
        end
    end
    
    %report - unseen
    for i = 9
        for j = 3 % "maintain"
             for k = 2
                 A{3}(:,:,i,j,k) = [0 0 0;%null
                                    1 1 0;%incorrect
                                    0 0 1];%correct
             end 
        end 
    end

    %report - seen 
    for i = 9
        for j = 3 % "maintain"
            for k = 3
                 A{3}(:,:,i,j,k) = [0 0 0;%null
                                    0 0 1;%incorrect
                                    1 1 0];%correct
            end 
        end 
    end 
    
    %--- Verbal report
    
    for i = 1:9
        for j = 1:3
            for k = 1:3
                    A{4}(:,:,i,j,k) = [1 1 1;%silent
                                       0 0 0;%"I"
                                       0 0 0;%"see"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end 
    
    % Report seen
    
    for i = 6
        for j = 3
            for k = 3
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       1 1 1;%"I"
                                       0 0 0;%"see"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 7
        for j = 3
            for k = 3
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       1 1 1;%"see"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 8
        for j = 3
            for k = 3
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       0 0 0;%"see"
                                       1 1 1;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 9
        for j = 3
            for k = 3
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       0 0 0;%"see"
                                       0 0 0;% "a"
                                       1 1 1;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    % Report not seen
    
    for i = 6
        for j = 3
            for k = 2
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       1 1 1;%"I"
                                       0 0 0;%"saw"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 7
        for j = 3
            for k = 2
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       0 0 0;%"saw"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       1 1 1;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 8
        for j = 3
            for k = 2
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       1 1 1;%"saw"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       0 0 0];%"anything"
            end 
        end 
    end
    
    for i = 9
        for j = 3
            for k = 2
                    A{4}(:,:,i,j,k) = [0 0 0;%silent
                                       0 0 0;%"I"
                                       0 0 0;%"saw"
                                       0 0 0;% "a"
                                       0 0 0;%"Gabor"
                                       0 0 0;%"didn't"
                                       1 1 1];%"anything"
            end 
        end 
    end

    a2 = A;

    %---- reset generative process (i.e. blanks are presented)

    for i = 4:9
          for j = 1:3
              for k = 1:3
                    A{1}(:,:,i,j,k) = [0 0 0;%left
                                       0 0 0;%right
                                       1 1 1];%blank
              end 
          end
    end


    % modulate precision of second level A-matrix to model WM gating
    for i = 2:3
            for k = 1:3
                %Low precision gating contents when policy is null & ~maintain 
                a2{1}(:,:,i,1:2,k) = spm_softmax(.5*log(a2{1}(:,:,i,1:2,k)+exp(-4)));
                a2{1}(:,:,i,3,k) = spm_softmax(1*log(a2{1}(:,:,i,3,k)+exp(-4)));   
            end 
    end

    %seperate generative process from generative model (multiplying by 64
    %prevents learning by making concentration parameteres very high).
    a2{1}= a2{1}*64;
    a2{2} = a2{2}*64;
    a2{3} = a2{3}*64;
    a2{4} = a2{4}*64;

    %transitions: B
    %--------------------------------------------------------------------------

    %Precision of sequence mapping
    B{1}(:,:,1) = eye(3,3); %Imprecise
    B{1}(:,:,2) = eye(3,3); %Precise

    B{2} = [0 0 0 0 0 0 0 0 0;
            1 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 1 1];


    %Working memory maintance
    B{3}(:,:,1) = [1 1 1;
                   0 0 0; 
                   0 0 0];
    B{3}(:,:,2) = [0 0 0
                   1 1 1;
                   0 0 0];   
    B{3}(:,:,3) = [0 0 0
                   0 0 0
                   1 1 1];   

    %Report
    B{4}(:,:,1) = [1 1 1;
                   0 0 0;
                   0 0 0];     
    B{4}(:,:,2) = [0 0 0;
                   1 1 1;
                   0 0 0];         
    B{4}(:,:,3) = [0 0 0;
                   0 0 0;
                   1 1 1];           

    b = B;

    %make b{1} of stimulus seqence hidden state factor imprecise
    b{1}(:,:,1) = spm_softmax(.5*log(B{1}(:,:,1)+exp(-4))); 
    b{1}(:,:,2) = spm_softmax(1.5*log(B{1}(:,:,2)+exp(-4))); 
    
    b{1} =  b{1}*100;
    b{2} =  b{2}*100;
    b{3} =  b{3}*100;
    b{4} =  b{4}*100;
     
    % Policies
    %--------------------------------------------------------------------------

     T = 9; %number of timesteps
     Nf = 4; %number of factors
     Pi = 3; %number of policies
     V = ones(T-1,Pi,Nf);

     %precise vs imprecise B-matrix
     V(:,:,1) = [2 2 2;
                 1 2 2;
                 1 2 2;
                 1 2 2;
                 1 2 2;
                 1 2 2;
                 1 2 2;
                 1 2 2];
     % ~maintain vs maintain
     V(:,:,3) = [1 1 1;
                 2 3 3;
                 2 3 3;
                 2 3 3;
                 2 3 3;
                 2 3 3;
                 2 3 3;
                 2 3 3];
    % Report: null = 1, unseen = 2, seen = 3
     V(:,:,4) = [1 1 1; 
                 1 1 1;
                 1 1 1;
                 1 1 1;
                 1 2 3;
                 1 2 3;
                 1 2 3;
                 1 2 3];
                    

    % C matrices (outcome modality by timesteps)
    %--------------------------------------------------------------------------
    C{1} = zeros(3,T);
    

    md = 1.5; % metabolic demand (i.e. how costly is working memory - higher = more costly)
    %metabolic demand of working memory
    C{2} = [0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 -md -md -md -md -md 0];
             
    if report == 1
        goal = 5;
    elseif report == 0
        goal = 0;
    end     
    
    %report
    C{3} = [0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 -1;
            0 0 0 0 0 0 0 0 goal];
        
    C{4} = zeros(7,T);

    % MDP Structure
    %--------------------------------------------------------------------------
    mdp.MDP  = MDP;
    mdp.link = [1 0 0 0; %lower level states by higher level observations
                0 1 0 0;
                0 0 0 1];

    mdp.T = T;                      % number of updates
    mdp.A = A;                      % observation model
    mdp.B = B;                      % transition probabilities
    mdp.C = C;                      % preferred outcomes
    mdp.D = D;                      % prior over initial states
    mdp.V = V;                      % policies

    mdp.d = d;
    mdp.a = a2;
    mdp.b = b;

    mdp.s = 1;                      % initial state
    mdp.erp = 1;
    mdp.alpha = 8;

    mdp.Aname = {'stimulus', 'Metabolic demand', 'Report Feedback', 'Verbal Report'};
    mdp.Bname = {'Sequence', 'Time in trial', 'WM goal', 'Report'};

    if cons == 1
        mdp_cons = spm_MDP_check(mdp);
    elseif cons == 2
        mdp_uncons = spm_MDP_check(mdp);
    end 

    clear A B C D V
    clear mdp
    clear MDP

end 
%% Construct MDP Structures
%==========================================================================

%conscious
rng('shuffle')
MDP_cons = spm_MDP_VB_X_PGNW(mdp_cons);

%unconscious
rng('shuffle')
MDP_uncons = spm_MDP_VB_X_PGNW(mdp_uncons);

%% MDP Plots
%==========================================================================

%conscious
spm_figure('GetWin','Seen - trial'); clf
spm_MDP_VB_trial(MDP_cons,4,4);
spm_figure('GetWin','Seen - ERP'); clf
spm_MDP_VB_ERP_PGNW(MDP_cons);
spm_figure('GetWin','Seen - LFP'); clf
spm_MDP_VB_LFP(MDP_cons,[],1);

%unconscious
spm_figure('GetWin','Unseen - trial'); clf
spm_MDP_VB_trial(MDP_uncons,4,4);
spm_figure('GetWin','Unseen - ERP'); clf
spm_MDP_VB_ERP_PGNW(MDP_uncons);
spm_figure('GetWin','Unseen - LFP'); clf
spm_MDP_VB_LFP(MDP_uncons,[],1);

%% Plot ERPs
%==========================================================================
[u1,~,ind] = spm_MDP_VB_ERP(MDP_cons,1); %Conscious
[u2,~,ind] = spm_MDP_VB_ERP(MDP_uncons,1); %Unconscious


i   = cumsum(ind);
i   = i(4) + (-127:-1);

u1  = u1(i,:);
u2  = u2(i,:);


pst = (1:length(u1));

limits = [1 86 -.3 .35];

figure(10)
hold on
plot(pst,sum(u1,2) ,'Color',[.7,0, 0], 'LineWidth',4) 
plot(pst,sum(u2,2) , 'Color',[.4,.4,.4], 'LineWidth',4) 
axis(limits)
ax = gca;
ax.YDir = 'reverse';
set(gca,'FontSize',30)
