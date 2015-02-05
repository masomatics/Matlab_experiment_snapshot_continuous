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
%%
    for(frame_index = 2:length(timesample))
        subplot(2, (length(timesample)-1)/2,frame_index-1)
        %scatter(snapshots(:,frame_index), hoge, 50,  tilde_pys(:,frame_index), 's');
    scatter(snapshots(1,:,frame_index),snapshots(2,:,frame_index))
        title(['t=', num2str(timesample(frame_index))], 'FontSize' ,20); 
    
        xlabel('x_1', 'FontSize', 20)
        ylabel('x_2', 'FontSize', 20)
    end