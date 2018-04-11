% **Basic models code base**           Jason Osik, 2012
% Following is a competitive decision-making network model between 
% two interconnected groups of neurons.  Initially the weights on connections 
% are set up with a combination of strong cross-inhibition and weak cross-excitation.
% Connections are modified between excitatory connections according to an STDP 
% learning rule and all connections are further modified by a synaptic
% scaling rule to maintain target firing rates.  Simulations are set to run
% over 8 trials, with the plasticity rules turned off for the first two
% trials to demonstrate the effects from reversing the sign of the input
% bias.  The input bias is again inverted during the last trial with
% plasticity in effect.

clear

dt = 0.0001;            %define time step
tmax = 5.0;             %set the full duration of each trial
t = 0:dt:tmax;          %create time vector

Ncells = 2;             %set the number of cellular assemblies
Ntrials = 8;            %set the number of simulation trials

Wrecurrent = 0.025;     %initial strength of recurrent excitatory connections
Wasym = -0.02;          %initial strength of cross-connections with strong net inhibition 
dW0 = 0.000005;             %initial magnitude of change on connection strengths

rmax0 = 1000;           %initial maximum firing rate
Ith0 = 8;               %current threshold for sigmoidal f-I activation function
Iwidth0 = 3;            %width parameter of sigmoidal function

tau_m = 0.010;          %time constant of membrane potential
tau_stdp = 0.020;       %time constant of stdp effect

Iapp0 = 0.2;            %strength of applied current due to stimulus
Idiff = 0.04;           %current bias 
Istart = 2.0;           %time of stimulus onset in seconds
Iend = 4.0;             %time that stimulus is turned off

epsilon = 0.00001;     %weight update parameter for synaptic scaling

W = zeros(Ncells,Ncells);   %initialize weight matrix
for cell1 = 1:Ncells,
    for cell2 = 1:Ncells,   %set all non-recurrent connections 
        W(cell1,cell2) = Wasym;
    end
    W(cell1,cell1) = Wrecurrent; %set all recurrent connections
end

for trial = 1:Ntrials,  %reverse bias between first and second trials then return it to original value; reverse bias again for final trial
    if (trial == 2)||(trial == 3)||(trial == Ntrials),
        Idiff = -Idiff;
    end
    
    r = zeros(length(t),Ncells);   %initialize rate matrix for full time course
    Iapp = zeros(length(t),Ncells); %initialize current input matrix
    
    imin = min(round(Istart/dt)+1,length(t)); %translate stimulus onset into time steps
    imax = min(round(Iend/dt)+1,length(t)); %translate stimulus offset into time steps
    
    for i = imin:imax,  %apply stimulus current over that duration
        for cell = 1:Ncells,
            Iapp(i,cell) = Iapp0 + Idiff*(cell-1)/Ncells;
        end
    end
    
    rmax = rmax0*ones(1,Ncells);
    Ith = Ith0*ones(1,Ncells);
    Iwidth = Iwidth0*ones(1,Ncells);
    
    for cell = 1:Ncells,
        r(1,cell) = 0.0; %set initial rate values to zero in all cells
    end
    
    for i = 2:length(t),    %integrate rate changes over full time course
        I = r(i-1,:)*W+Iapp(i,:);
        rinf = rmax./(1.+exp(-(I-Ith)./Iwidth));
        r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau_m);
        
    end
    
    spikes = zeros(length(t),Ncells); %initialize spike matrix
    for cell = 1:Ncells,    %generate spikes in each cell according to Poisson process 
        for i = 1:length(t),
            if r(i,cell)*dt > rand(1),
                spikes(i,cell) = 1;
            end
        end
    end
    
    figure()    %plot rates in each cell
    for cell = 1:Ncells,
        subplot(Ncells,1,cell);
        plot(t,r(:,cell));
        xlabel('Time, in secs.');
        ylabel('Rate, in Hz');
        axis([0 tmax 0 rmax0]);
    end         %display weight matrix
        figure()
        imagesc(W);
        xlabel('Presynaptic');
        ylabel('Postsynaptic');
        colorbar;
        
    % the following high-level loop contains plasticity and homeostatic elements
    if (trial > 2),         %do not apply plasticity to first 2 trials, to demonstrate effect of bias alone
        for cell1 = 1:Ncells,       %for all presynaptic cells
        
            st1 = find(spikes(:,cell1));    %find spike times
            N1 = length(st1);               %total number of spikes
            for cell2 = 1:Ncells,   %for all postsynaptic cells (each pre-post pair defines a connection)
                dwsum = 0.0;        %start with summed weight change for the given connection set to zero
            
                if cell1 ~= cell2,  %update only off-diagonal (i.e., non-recurrent) connections 
                    st2 = find(spikes(:,cell2));    %find postsynaptic spike times
                    N2 = length(st2);       %total number of postsynaptic spikes
                    for i = 1:N1    %integrating over entire presynaptic spike train
                        for j = 1:N2 %in conjunction with postsynaptic spike train
                            deltat = (st2(j)-st1(i))*dt;    %compute intervals between pre/post spikes
                            if (deltat > 0),    %if pre-before-post
                                dwsum = dwsum + exp(-deltat/tau_stdp);  %increase weight
                            end
                            if (deltat < 0),    %if post-before-pre
                                dwsum = dwsum - exp(deltat/tau_stdp);   %decrease weight
                            end
                        end
                    end
                    dwsum = dwsum*dW0,  %scale weight change
                    dwsum = min(dwsum,0.5); %and set upper and lower bounds
                    dwsum = max(dwsum,-0.5);
                    W(cell1,cell2) = W(cell1,cell2)*(1 + dwsum);%complete the update to the appropriate element of the weight matrix
                    W,
                end
            end
        end
            %compute active versus non-active (baseline) firing rates
            %(viewing every cell in terms of postsynaptic activity)
        for k = 1:Ncells,       %determine the baseline postsynaptic firing rate before stimulus presentation
            r_pre_mean(1:k) = mean(r(1,k):r(imin,k));  %mean firing rate between t=0 and time of stimulus
        end
        r_goal = r_pre_mean;    %make this the postsynaptic goal mean rate
        
        for k = 1:Ncells,       %determine the mean postsynaptic firing rate during stimulus presentation
            r_stim_mean(1:k) = mean(r(imin,k):r(imax,k)); %mean firing rate during duration of stimulus
        end
        r_act = r_stim_mean;    %set the active mean postsynaptic firing rate equal to this
        
        for i = 1:Ncells,   %update weights according to multiplicative synaptic scaling
            for j = 1:Ncells,
                if (W(i,j) > 0),    %implement homeostasis on net excitatory elements
                    W(i,j) = W(i,j)*(1 - epsilon*(r_act(1,j) - r_goal(1,j)));
                end
                if (W(i,j) < 0),    %implement homeostasis on net inhibitory elements
                    W(i,j) = W(i,j)*(1 + epsilon*(r_act(1,j) - r_goal(1,j)));
                end
            end
        end
        
    
    end
    
end

%%
% Alternatively, the homeostasis block above can be written:
%        for i = 1:Ncells,   %update weights according to multiplicative synaptic scaling
%            for j = 1:Ncells,
%                delta_r = (r_act(1,j) - r_goal(1,j));
%                if (delta_r > 0),    %implement homeostasis on net excitatory elements
%                    W(i,j) = W(i,j)*(1 - epsilon*(delta_r(1,j));
%                end
%                if (delta_r < 0),    %implement homeostasis on net inhibitory elements
%                    W(i,j) = W(i,j)*(1 + epsilon*(delta_r(1,j));
%                end
%            end
%        end
%   

