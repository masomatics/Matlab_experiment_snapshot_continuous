theta(1) = 2 ; 
theta(2) = 10;
theta(3) = 1/4;
theta(4) = 1;
tend = 5; 

sigV = 2;sigW = 2;
initx = [1; 1]

num_timepts = 25;
Ntry = 3;

rnsource = randn([2, Ntry, num_timepts]);

[timepts,datapts, derivpts, derivpts2] = and_CFD_datagen_mass_deriv2(initx, tend, theta, sigV, num_timepts, rnsource, Ntry);

timesample = 1:1:5;
timeindex = find(ismember(timepts, timesample));

rnsource2 = randn([2, Ntry, length(timesample)]);
snapshots = datapts(:, :, timeindex) + rnsource2;

theta0= rand(1,4).^2;
%%
rnsource = randn([2, Ntry, num_timepts]);

[datmat, tilde_pys, deriv1] = and_CFD_datagen_mass_derivStat_all_parameters...
                        (initx, tend, theta0, sigV, sigW, num_timepts, ...
                        rnsource, snapshots, timesample, Ntry);
                
[datmat, tilde_pys, deriv2] = and_CFD_datagen_mass_derivStat_all_parameters_temp...
                        (initx, tend, theta0, sigV, sigW, num_timepts,...
                        rnsource, snapshots, timesample, Ntry)  ;                 
deriv1 
deriv2
%%
close all;
N = 2
rnsource_debug =  randn([2, N, num_timepts]);
%%
close all;
[datmat, tilde_pys, deriv] = and_CFD_datagen_mass_derivStat_all_parameters_temp...
                        (initx, tend, theta0, sigV, sigW, num_timepts,...
                        rnsource_debug, snapshots, timesample, N);

close all;
%%
[datmat, tilde_pys, deriv1] = and_CFD_datagen_mass_derivStat_all_parameters...
                        (initx, tend, theta0, sigV, sigW, num_timepts, ...
                        rnsource_debug, snapshots, timesample, N);

close all;