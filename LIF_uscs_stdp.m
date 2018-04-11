%  **Basic models code base**           Jason Osik, 2012
%  LIF neuron model w/ 2 excitatory synaptic inputs representing
%  conditioned and unconditioned stimuli where synaptic strengths are
%  updated according to a spike-timing dependent plasticity rule.  

clear
clf

Vth = -0.054;
V_L = -0.070;
Vreset = -0.080;
Vex = 0.0;
c_m = 1e-8;
g_L = 1e-6;
tau_m = c_m/g_L;

Ninputs = 2;
Ntrials = 14;   % Total trials; the sum of the trial subset sizes below must equal this value
Ntrialset1 = 5; % Size of trial subset 1
Ntrialset2 = 2; % Size of trial subset 2
Ntrialset3 = 5; % Size of trial subset 3
Ntrialset4 = 2; % Size of trial subset 4
TrialNames1 = (1:Ntrialset1);
TrialNames2 = (Ntrialset1+1:Ntrialset1+Ntrialset2);
TrialNames3 = (TrialNames2(end)+1:TrialNames2(end)+Ntrialset3);
TrialNames4 = (TrialNames3(end)+1:TrialNames3(end)+Ntrialset4);

tmax = 0.200;   % time for one trial run
dt = 0.001;
t = 0:dt:tmax;

sigma = 1.0e-7;

tau_syn = 0.005;
A_LTP = 0.35e-6;
A_LTD = 0.45e-6;
tau_ltp = 0.025;
tau_ltd = 0.035;
g_US = 1.2e-6;
g_CS = 0.0;
g_CSmax = 1.2e-6;


noise = randn(1,length(t))*sigma*dt;

for trial = 1:Ntrials,
    s1 = zeros(size(t));
    s2 = zeros(size(t));
    
    if (trial >= TrialNames1(1))&&(trial <= TrialNames1(end)),
        US_onset = 0.100;
        CS_onset = 0.090;
    end
    if (trial >= TrialNames2(1))&&(trial <= TrialNames2(end)),
        US_onset = 2*tmax*Ntrials;
        CS_onset = 0.090;
    end
    if (trial >= TrialNames3(1))&&(trial <= TrialNames3(end)),
        US_onset = 0.100;
        CS_onset = 0.110;
    end
    if (trial >= TrialNames4(1))&&(trial <= TrialNames4(end)),
        US_onset = 2*tmax*Ntrials;
        CS_onset = 0.110;
    end

   I_tot = zeros(length(t)); 
   spikesin = zeros(length(t),Ninputs);
   spikesout = zeros(size(t));
   V = zeros(size(t));
   V(1) = V_L;
    
   for i = 2:length(t),
        if (i-1 == CS_onset/dt), 
            spikesin(i-1,2) = 1;
            s2(i-1) = s2(i-1) + 1;
            
        end
        if (i-1 == US_onset/dt),  
            spikesin(i-1,1) = 1;
            s1(i-1) = s1(i-1) + 1;
        
        end
        
        I_tot(i-1) = I_tot(i-1) - (g_US*s1(i-1) + g_CS*s2(i-1))*(V(i-1)-Vex);
        V_inf = V_L + (I_tot(i-1))/g_L;
        V(i) = V_inf + (V(i-1)-V_inf)*exp(-dt/tau_m);
        
        s1(i) = s1(i-1)*exp(-dt/tau_syn);
        s2(i) = s2(i-1)*exp(-dt/tau_syn);
        
        if V(i) >= Vth,
           spikesout(i) = 1;
           V(i) = Vreset;
        end
    end
   
    figure(trial)
    plot(t,V);
    xlabel('t, in secs.');
    ylabel('Post-synaptic membrane potential (in V)');
    drawnow
    
    series_out = find(spikesout);
    N_output = length(series_out);
    series_CSin = find(spikesin(:,2));
    N_CSin = length(series_CSin);
     
    for i = 1:N_output,
        for j = 1:N_CSin,
            delta_t = dt*(series_out(i)-series_CSin(j));
            if (delta_t > 0)
                g_CS = g_CS + A_LTP*exp(-delta_t/tau_ltp);
            end
            if (delta_t < 0)
                g_CS = g_CS - A_LTD*exp(-delta_t/tau_ltd);
            end
        end
    end
    
    g_CS = min(g_CS,g_CSmax);
    g_CS = max(g_CS,0);
    
    g_CSvar(trial) = g_CS;

end
  
x = 1:trial;
figure;
plot(x,g_CSvar,'k-');
xlabel('Trial number');
ylabel('CS synaptic conductance, g_CS');
drawnow