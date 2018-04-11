%  **Basic models code base**           Jason Osik, 2012
%Behavior of a Rate model for orientation tuning based on
%Ben-Yishai, Bar-Or, & Sompolinsky, PNAS (USA) 92:3844 (1995)
%%
%1.

dt = 0.005;             %set time step
tmax = 5.0;             %set duration of simulation
t = 0:dt:tmax;          %create time vector

Ncells = 50;            %define the number of cells in the simulation
dtheta = pi/Ncells;     %define size of graduations in orientation angle
r = zeros(length(t),Ncells); %initialize rate matrix for all cells at all time points
r(1,:) = 10*rand(1,Ncells);  %start firing rate of all cells at random value between 0-10Hz
h = zeros(1,Ncells);    %initialize thalamic input vector
I_total = zeros(1,Ncells); %initialize the total input vector
rmax = 300;             %set max. firing rate **NOTE: this is in effect but not generally
                        %used in our simulations at the parameter values
                        %given

A = 40.;                %set max. thalamic input rate
epsilon = 0.1;          %set cortical anisotropy (related to aspect ratio) of input
J0 = -0.5;              %set net inhibitory component
J2 = 1;                 %set net orientation-specific component
tau = 0.050;            %set time constant
c = [0.1;0.2;0.4;0.8];  %define vector of stimulus contrast levels to be tested
cueCell = 25;   % Ncells vector index for cell with theta = 0 (cue theta)

for x = 1:length(c),    %test over all contrast levels
    for cell = 1:Ncells,    %define the total thalamic input to all cells
        h(cell) = A*c(x)*(1 - epsilon + epsilon*cos(2*pi*(cell-cueCell)/Ncells));
    end

    for i = 2:length(t),    %integrate over all time steps
        for cell_j = 1:Ncells,  % define feedforward input component to all cells
            I_total(cell_j) = h(cell_j);
            for cell_k = 1:Ncells,  % define total output to all other cells as combination of feedforward and recurrent input
                I_total(cell_j) = I_total(cell_j)+r(i-1,cell_k)*(1/pi)*dtheta*(J0+J2*cos(2*pi*(cell_j-cell_k)/Ncells));
            end
        end
    

        for cell = 1:Ncells,    % set the steady-state firing rates for all cells according to threshold linear f-I curve (w/ hard upper bound)
            if I_total(cell) > 0&&I_total(cell) < rmax,
                rinf(cell) = I_total(cell);
            elseif I_total(cell) > rmax,
                rinf(cell) = rmax;
            elseif I_total(cell) <= 0,
                rinf(cell) = 0;
            end
        end
            % evaluate firing rate for all cells at next time step
    r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau);
    end
    %translate cell numbers to preferred orientations (first in radians
    %then degrees
theta = -(pi/2)+dtheta:dtheta:pi/2;
theta_deg = 57.2958*theta;
    % plot activity profile at all values of c
figure(1);
    if x == 1,
        plot(theta_deg,r(end,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,r(end,:),'r');
    elseif x == 3,
        plot(theta_deg,r(end,:),'g');
    elseif x == 4,
        plot(theta_deg,r(end,:),'b');
    end
    % plot thalamic input profile at all values of c   
figure(2);
    if x == 1,
        plot(theta_deg,h(1,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,h(1,:),'r');
    elseif x == 3,
        plot(theta_deg,h(1,:),'g');
    elseif x == 4,
        plot(theta_deg,h(1,:),'b');
    end

end
% label appropriate figures
figure(1);
hold on;
xlabel('Orientation (in degrees)');
ylabel('Firing rate (Hz)');
legend('c=0.1','c=0.2','c=0.4','c=0.8');
hold off;
figure(2);
hold on;
xlabel('Orientation (in degrees)');
ylabel('h, Thalamic input');
legend('c=0.1','c=0.2','c=0.4','c=0.8');



%%
%2.  The same model is examined with the following parameter modifications:
% J0 = -7.3, J2 = 11, A = 40 Hz, epsilon = 0.1
%NOTE: %*** - indicates lines that have been changed from part 1.
 
dt = 0.005;             
tmax = 5.0;             
t = 0:dt:tmax;          

Ncells = 50;            
dtheta = pi/Ncells;     
r = zeros(length(t),Ncells); 
r(1,:) = 10*rand(1,Ncells);  
h = zeros(1,Ncells);    
I_total = zeros(1,Ncells);
rmax = 300;

A = 40.;
epsilon = 0.1;
J0 = -7.3;      %***
J2 = 11;        %***
tau = 0.050;
c = [0.1;0.2;0.4;0.8];
cueCell = 25;   

for x = 1:length(c),
    for cell = 1:Ncells,
        h(cell) = A*c(x)*(1 - epsilon + epsilon*cos(2*pi*(cell-cueCell)/Ncells));
    end

    for i = 2:length(t),
        for cell_j = 1:Ncells,
            I_total(cell_j) = h(cell_j);
            for cell_k = 1:Ncells,
                I_total(cell_j) = I_total(cell_j)+r(i-1,cell_k)*(1/pi)*dtheta*(J0+J2*cos(2*pi*(cell_j-cell_k)/Ncells));
            end
        end
    

        for cell = 1:Ncells,
            if I_total(cell) > 0&&I_total(cell) < rmax,
                rinf(cell) = I_total(cell);
            elseif I_total(cell) > rmax,
                rinf(cell) = rmax;
            elseif I_total(cell) <= 0,
                rinf(cell) = 0;
            end
        end

    r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau);
    end

theta = -(pi/2)+dtheta:dtheta:pi/2;
theta_deg = 57.2958*theta;

figure(3);
    if x == 1,
        plot(theta_deg,r(end,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,r(end,:),'r');
    elseif x == 3,
        plot(theta_deg,r(end,:),'g');
    elseif x == 4,
        plot(theta_deg,r(end,:),'b');
    end
    
figure(4);
    if x == 1,
        plot(theta_deg,h(1,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,h(1,:),'r');
    elseif x == 3,
        plot(theta_deg,h(1,:),'g');
    elseif x == 4,
        plot(theta_deg,h(1,:),'b');
    end

end

figure(3);
hold on;
xlabel('Orientation (in degrees)');
ylabel('Firing rate (Hz)');
legend('c=0.1','c=0.2','c=0.4','c=0.8');
hold off;
figure(4);
hold on;
xlabel('Orientation (in degrees)');
ylabel('h, Thalamic input');
legend('c=0.1','c=0.2','c=0.4','c=0.8');

%%
%3. Using the same model, we consider how randomness in the distribution of
% thalamic input, h, influences activity
%Part 1 

dt = 0.005;            
tmax = 5.0;             
t = 0:dt:tmax;          

Ncells = 50;            
dtheta = pi/Ncells;     
r = zeros(length(t),Ncells); 
r(1,:) = 10*rand(1,Ncells);  
h = zeros(1,Ncells);    
I_total = zeros(1,Ncells);
rmax = 300;
sigma_h = 1.;
mean_eta = 0;
sd_eta = 1;


A = 40.;
epsilon = 0.1;
J0 = -0.5;  %***
J2 = 1;     %***
tau = 0.050;
c = [0.1;0.2;0.4;0.8];
cueCell = 25;   

for x = 1:length(c),
    
eta_i = mean_eta+(sd_eta*randn(1,Ncells));  %*** define Gaussian component of noise term
    
    for cell = 1:Ncells,
        h(cell) = A*c(x)*(1 - epsilon + epsilon*cos(2*pi*(cell-cueCell)/Ncells))+...
            (sigma_h*eta_i(cell));      %*** add total noise term to thalamic input
    end

    for i = 2:length(t),
        for cell_j = 1:Ncells,
            I_total(cell_j) = h(cell_j);
            for cell_k = 1:Ncells,
                I_total(cell_j) = I_total(cell_j)+r(i-1,cell_k)*(1/pi)*dtheta*(J0+J2*cos(2*pi*(cell_j-cell_k)/Ncells));
            end
        end
    

        for cell = 1:Ncells,
            if I_total(cell) > 0&&I_total(cell) < rmax,
                rinf(cell) = I_total(cell);
            elseif I_total(cell) > rmax,
                rinf(cell) = rmax;
            elseif I_total(cell) <= 0,
                rinf(cell) = 0;
            end
        end

    r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau);
    end

theta = -(pi/2)+dtheta:dtheta:pi/2;
theta_deg = 57.2958*theta;

figure(5);
    if x == 1,
        plot(theta_deg,r(end,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,r(end,:),'r');
    elseif x == 3,
        plot(theta_deg,r(end,:),'g');
    elseif x == 4,
        plot(theta_deg,r(end,:),'b');
    end
    
figure(6);
    if x == 1,
        plot(theta_deg,h(1,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,h(1,:),'r');
    elseif x == 3,
        plot(theta_deg,h(1,:),'g');
    elseif x == 4,
        plot(theta_deg,h(1,:),'b');
    end

end

figure(5);
hold on;
xlabel('Orientation (in degrees)');
ylabel('Firing rate (Hz)');
legend('c=0.1','c=0.2','c=0.4','c=0.8');
hold off;
figure(6);
hold on;
xlabel('Orientation (in degrees)');
ylabel('h, Thalamic input');
legend('c=0.1','c=0.2','c=0.4','c=0.8');

%%
%3.
%Part 2 Using parameter changes made for 2.

dt = 0.005;            
tmax = 5.0;             
t = 0:dt:tmax;          

Ncells = 50;           
dtheta = pi/Ncells;     
r = zeros(length(t),Ncells); 
r(1,:) = 10*rand(1,Ncells); 
h = zeros(1,Ncells);    
I_total = zeros(1,Ncells);
rmax = 300;
sigma_h = 4.;
mean_eta = 0;
sd_eta = 1;


A = 40.;
epsilon = 0.1;
J0 = -7.3;  %***
J2 = 11;    %***
tau = 0.050;
c = [0.1;0.2;0.4;0.8];
cueCell = 25;   

for x = 1:length(c),
    
eta_i = mean_eta+(sd_eta*randn(1,Ncells));    %***
    
    for cell = 1:Ncells,
        h(cell) = A*c(x)*(1 - epsilon + epsilon*cos(2*pi*(cell-cueCell)/Ncells))+...
            (sigma_h*eta_i(cell));            %***
    end

    for i = 2:length(t),
        for cell_j = 1:Ncells,
            I_total(cell_j) = h(cell_j);
            for cell_k = 1:Ncells,
                I_total(cell_j) = I_total(cell_j)+r(i-1,cell_k)*(1/pi)*dtheta*(J0+J2*cos(2*pi*(cell_j-cell_k)/Ncells));
            end
        end
    

        for cell = 1:Ncells,
            if I_total(cell) > 0&&I_total(cell) < rmax,
                rinf(cell) = I_total(cell);
            elseif I_total(cell) > rmax,
                rinf(cell) = rmax;
            elseif I_total(cell) <= 0,
                rinf(cell) = 0;
            end
        end

    r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau);
    end

theta = -(pi/2)+dtheta:dtheta:pi/2;
theta_deg = 57.2958*theta;

figure(7);
    if x == 1,
        plot(theta_deg,r(end,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,r(end,:),'r');
    elseif x == 3,
        plot(theta_deg,r(end,:),'g');
    elseif x == 4,
        plot(theta_deg,r(end,:),'b');
    end
    
figure(8);
    if x == 1,
        plot(theta_deg,h(1,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,h(1,:),'r');
    elseif x == 3,
        plot(theta_deg,h(1,:),'g');
    elseif x == 4,
        plot(theta_deg,h(1,:),'b');
    end

end

figure(7);
hold on;
xlabel('Orientation (in degrees)');
ylabel('Firing rate (Hz)');
legend('c=0.1','c=0.2','c=0.4','c=0.8');
hold off;
figure(8);
hold on;
xlabel('Orientation (in degrees)');
ylabel('h, Thalamic input');
legend('c=0.1','c=0.2','c=0.4','c=0.8');

%%
%4.
% We consider the original model settings, but with no thalamic tuning
% (i.e., epsilon = 0) together with varying levels of J0 and J2, which
% correspond to different degrees of excitation-inhibition balance and
% different levels of recurrent signal amplitude, respectively (Ben-Yishai
% et al., pg. 3845.

% We examine J0 values of -1 in conjuction with J2 values of 1,1.5,2.5,4.5

dt = 0.005;             
tmax = 5.0;             
t = 0:dt:tmax;          

Ncells = 50;            
dtheta = pi/Ncells;     
r = zeros(length(t),Ncells); 
r(1,:) = 10*rand(1,Ncells);  
h = zeros(1,Ncells);    
I_total = zeros(1,Ncells);
rmax = 300;

A = 40.;
epsilon = 0.0;  %*** Effect of thalamic input on orientation (aspect ratio) is set to zero
J0 = -1;    %***
J2 = 4.5;   %***    These values were varied to multiple different combinations
tau = 0.050;
c = [0.1;0.2;0.4;0.8];
cueCell = 25;   

for x = 1:length(c),
    for cell = 1:Ncells,
        h(cell) = A*c(x)*(1 - epsilon + epsilon*cos(2*pi*(cell-cueCell)/Ncells));
    end

    for i = 2:length(t),
        for cell_j = 1:Ncells,
            I_total(cell_j) = h(cell_j);
            for cell_k = 1:Ncells,
                I_total(cell_j) = I_total(cell_j)+r(i-1,cell_k)*(1/pi)*dtheta*(J0+J2*cos(2*pi*(cell_j-cell_k)/Ncells));
            end
        end
    

        for cell = 1:Ncells,
            if I_total(cell) > 0&&I_total(cell) < rmax,
                rinf(cell) = I_total(cell);
            elseif I_total(cell) > rmax,
                rinf(cell) = rmax;
            elseif I_total(cell) <= 0,
                rinf(cell) = 0;
            end
        end

    r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tau);
    end

theta = -(pi/2)+dtheta:dtheta:pi/2;
theta_deg = 57.2958*theta;

figure(9);
    if x == 1,
        plot(theta_deg,r(end,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,r(end,:),'r');
    elseif x == 3,
        plot(theta_deg,r(end,:),'g');
    elseif x == 4,
        plot(theta_deg,r(end,:),'b');
    end
    
figure(10);
    if x == 1,
        plot(theta_deg,h(1,:),'k');
        hold on;
    elseif x == 2,
        plot(theta_deg,h(1,:),'r');
    elseif x == 3,
        plot(theta_deg,h(1,:),'g');
    elseif x == 4,
        plot(theta_deg,h(1,:),'b');
    end

end

figure(9);
hold on;
xlabel('Orientation (in degrees)');
ylabel('Firing rate (Hz)');
legend('c=0.1','c=0.2','c=0.4','c=0.8');
hold off;
figure(10);
hold on;
xlabel('Orientation (in degrees)');
ylabel('h, Thalamic input');
legend('c=0.1','c=0.2','c=0.4','c=0.8');