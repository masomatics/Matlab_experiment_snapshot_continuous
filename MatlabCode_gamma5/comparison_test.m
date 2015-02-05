hoge = load('snapshots_2parameters.mat');
snapshots = hoge.snapshots ;
theta0 = [1/4,10];

%%
close all;
Ntry = 5;
%tic,
rnsource = randn([2, Ntry, num_timepts]);
%toc
timesample = [     1     2     3     4     5];

smallsnap = snapshots(:, 1:6, :) 
snapshots = smallsnap;


num_timepts = 2500;
sigV = 2;
sigW = 2;
tend = 5;
initx = [1 ; 0];

%%
theta0_aug = [2, theta0(2), theta0(1), 1]
[datmat2, tilde_pys2, deriv2] = and_CFD_datagen_mass_derivStat_all_parameters...
                        (initx, tend, theta0_aug, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);
                    
%%                    
[datmat, tilde_pys, deriv] = and_CFD_datagen_mass_derivStat2...
                        (initx, tend, theta0, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);
                    
          datmat=          permute((datmat(1,1,:)),[2,3,1]);
          datmat2=         permute((datmat2(1,1,:)),[2,3,1]);
          
          %%
          figure(100)
          plot(datmat, 'b')
          hold on;
          plot(datmat2,'g')
          
          max((datmat - datmat2).^2)
          
          hold off;
          
          figure(101)
          hist(tilde_pys(2,:)',300) 
          hold on;
          hist(tilde_pys2(2,:)',300) 