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
num_iter = 1500;

%theta0 = [0.11054      7.1887];
%theta0 = [0.16072      8.2049]
%theta0 = [1.9474 9  0.2331 0.8698]; 

%theta0= rand(1,4).^2;
theta0 = [1.9578      8.3635     0.25871     0.74678]; %1600th
theta0 = [1.9075      8.6019     0.39438     0.92227]; %1750th

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


