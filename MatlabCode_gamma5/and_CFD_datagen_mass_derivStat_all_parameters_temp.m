
%
% This is the diffusion version of Dave's CFD paper. 
%
% dXt = A Xt dt + dW.  
% 
% rnsource  2 by  sample number by numtimepoints 
% dataset   2 by  sample number by numtimepoints 
% init      2 by 1
%
% snapshots snapshots.  Each column corresponds to a time.
% snaptime  sampled time. length is num_frames. 
%

function [datmat, tilde_pys, derivative] = and_CFD_datagen_mass_derivStat_all_parameters...
                        (init, tend, theta, sigV, sigW, num_timepts, rnsource, snapshots, snaptime, N)
                    
            include= [];
            
%%%preset variables %%%                     
    [species, num_particles, num_slices] = size(snapshots);  
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    num_parameters = 4;    
    rxn_mat =   [1 0 -1 0; 0 1 0 -1]; 
    %This part is specific for "ALL particles case"
    compress_snap_wgts = 1/num_particles * ones(num_particles, num_slices);    
    
%%%Target variables%%%     
    derivative = zeros(num_parameters,1);
    tilde_p_ymk = zeros(num_particles);
    p_ymkj = zeros(num_particles, N);
    dp_jr = zeros(N,num_particles);
    dEP_kr= zeros(N, num_particles); 
%%%Temporary variables%%%    
    data_now = NaN(species, N); 
    rxn_rate = zeros(4,N);
    ftheta = NaN(num_parameters, N); 
    ycopy = zeros(species, num_particles, N); 
    xcopy = zeros(species, num_particles, N);
    scale_p_ymkj = zeros(1,N) ;
    scale_p_ymkj_mat = zeros(num_particles,N);
    derivative_alltime = zeros(num_parameters, num_slices); 
 
%%%Initialization%%%
    v = sigV*rnsource;
    initx = repmat(init,1,N);
    data_now = initx;
    timemat = delta*(0:(num_timepts+1));    
    deriv_loglike = zeros(num_parameters,N);
    snaptime_now  = 1;

    
%%%Debug : Assessment variables %%%             
    energy = zeros(1, num_slices);  

%%%Monte Carlo Simulation%%%     
    for(m = 2 : (num_timepts+1))
        

        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * data_now(1,:);
        rxn_rate(3,:) = theta(3) * data_now(1,:);
        rxn_rate(4,:) = theta(4) * data_now(2,:);        
        
        waitbar(m/num_timepts);
        %%% ftheta part is considered here.
        deriv_loglike(1,:) = deriv_loglike(1,:) +  rxn_mat(:,1)'*(v(:,:,m-1)*sqrt(delta))/(sigV^2);
        deriv_loglike(2,:) = deriv_loglike(2,:) +  data_now(1,:).*(rxn_mat(:,2)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)   ;  
        deriv_loglike(3,:) = deriv_loglike(3,:) +  data_now(1,:).*(rxn_mat(:,3)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)  ;
        deriv_loglike(4,:) = deriv_loglike(4,:) +  data_now(2,:).*(rxn_mat(:,4)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)  ;        
        
        %Update the data
        data_now = data_now + rxn_mat*rxn_rate*delta + sqrt(delta)*v(:,:,m-1);

        
        if(timemat(m) == snaptime(snaptime_now))

            %%% New Codes %%% 
            ycopy = repmat(snapshots(:,:,snaptime_now), [1,1,N]);
            xcopy = permute(repmat(data_now,[1,1,num_particles])...
                ,[1,3,2]);
             distance_pair = min(squeeze(sum((xcopy - ycopy).^2,1)/sigW^2), 300);
            p_ymkj = exp(-(distance_pair));
            scale_p_ymkj = 1./sum(p_ymkj,1);
%            scale_p_ymkj_mat = diag(scale_p_ymkj); %Multiplying diag will
%            cost O(N^3). Elementwise will cost O(N^2).  
            scale_p_ymkj_mat = repmat(scale_p_ymkj,[num_particles,1]) ;
    %tilde_pym in the old version
            p_ymkj = p_ymkj .* scale_p_ymkj_mat;
            tilde_p_ymk = mean(p_ymkj,2);
            dEP_kr = 1/N *p_ymkj * deriv_loglike' ;
            
            derivative_alltime(:,snaptime_now) = compress_snap_wgts(:,snaptime_now)' * ...
                (dEP_kr ./ repmat(tilde_p_ymk,[1, num_parameters]));

            snaptime_now = snaptime_now  +1;

                   if(max(p_ymkj(:))  > min(p_ymkj(:)) )
                        %In the case of using empirical distribution,
                        %the true distribution can be soo far from the
                        %simulated distribution that the energy can be
                        %wrongfully too high (all uniform) This part is
                        %very hand-wavy.... I am ridding of such frames
                        %from the consideration.                         
                        include = [include, snaptime_now];
                   end
        end %end of frame
                               

    end %end of time loop 
    derivative = mean(derivative_alltime(:,(include-1)),2);
    datmat = 0;
    tilde_pys = 0;
    close(h);
end 

function [rxn_rate] = compute_rxn_rate(theta, data)
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * data(1,:);
        rxn_rate(3,:) = theta(3) * data(1,:);
        rxn_rate(4,:) = theta(4) * data(2,:);
end 
