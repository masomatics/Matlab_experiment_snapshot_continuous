%%
clear all;
theta(1) = 2 ; 
theta(2) = 10;
theta(3) = 1/4;
theta(4) = 1;
tend = 5; 

sigV = 2;
sigW = 2;
initx = [0; 0]
%num_timepts = 2500;
num_timepts = 2500;
%Ntry = 10000;
Ntry = 5000;

rnsource = randn([2, Ntry, num_timepts]);

[timepts,datapts, derivpts, derivpts2] = and_CFD_datagen_mass_deriv2(initx, tend, theta, sigV, num_timepts, rnsource, Ntry);

timesample = 1:1:5;
timeindex = find(ismember(timepts, timesample));

rnsource2 = randn([2, Ntry, length(timesample)]);
snapshots = datapts(:, :, timeindex) + rnsource2;

%save('snapshots_allparameters.mat', 'snapshots')

energyhistory_cum + []
%%
hoge = load('snapshots_allparameters.mat');
snapshots = hoge.snapshots;

close all;
Ntry = 5000;
%tic,
rnsource = randn([2, Ntry, num_timepts]);
%toc

eta = 0.05; 
%num_iter = 1600;
num_iter = 500;

%theta0 = [0.11054      7.1887];
%theta0 = [0.16072      8.2049]
%theta0 = [1.9474 9  0.2331 0.8698]; 
%timesample = [ 1     2     3     4     5];
%theta0= rand(1,4).^2;
theta0 = [1.9578      8.3635     0.25871     0.74678]; %1600th
theta0 = [1.9075      8.6019     0.39438     0.92227]; %1750th
theta0 = [1.9664      8.9375     0.23683     0.84142]; %3250th
theta0 = [1.9648      8.9795     0.23459     0.84587]; %3550th
theta0 = [1.9704      9.0167     0.24467     0.84785]; %3900th
theta0 = [1.9628       9.032     0.23819     0.85425]; %4100th
theta0 = [1.9696      9.0521     0.24324     0.83667]; %4300th
theta0 = [1.9626      9.0638     0.24049     0.85051]; %4500th
theta0 = [1.959      9.1515      0.2433      0.86486]; %6500th
theta0 = [1.9582      9.1592     0.24799     0.87206]; %7000th

energyhistory = zeros(1, num_iter) ;


thetahistory = zeros(1+num_iter,length(theta0));
thetahistory(1,:) = theta0;


for iter = 1:num_iter
    rnsource = randn([2, Ntry, num_timepts]);

    display(['initiating the ', num2str(iter), 'th iterations...']);
    tic,
    [derivative, energy] = and_CFD_datagen_mass_derivStat_all_parameters_totVar...
                        (initx, tend, theta0, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);
    display([num2str(iter), 'th iteration : deriv= ' , num2str(derivative')]);
    

    
    energyhistory(iter) = sum(energy); 

    
    %I am trying to increase the energy, so I will take positive gradient
    theta0 = theta0 - eta*derivative';
    theta0 = max(theta0, 0);
    toc
        display([num2str(iter), 'th iteration complete:  new theta0 is =' , num2str(theta0)]);
        display([num2str(iter), 'th iteration complete:  new energy is =' , num2str(sum(energy))]);

        thetahistory(iter+1,:) = theta0;
end 

energyhistory_cum = [energyhistory_cum, energyhistory];

