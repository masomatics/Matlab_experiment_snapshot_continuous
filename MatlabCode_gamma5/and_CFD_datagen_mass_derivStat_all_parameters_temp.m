
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
    
    [species, num_particles, num_frames] = size(snapshots);  
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    num_parameter = 4;    
    rxn_mat =   [1 0 -1 0; 0 1 0 -1]; 
    rxn_rate = zeros(4,N);    
    timemat = delta*(0:(num_timepts+1));    
    datmat = zeros(2, N, num_timepts +1);    
    %\beta_m^k
    derivdat = NaN(num_particles,num_frames);    
    deriv_loglike = zeros(num_parameter,N);
   % deriv_loglike2 = zeros(5,N);    
    deriv_pym = zeros(num_parameter, N, num_particles);
    tilde_pym = zeros(N, num_particles);
    sum_tilde_pym = zeros(N,1);    
    deriv_py = zeros(num_parameter,num_particles);
    tilde_py = zeros(1,num_particles);
    tilde_pys= zeros(num_frames, num_particles);    
    deriv_pm = zeros(num_parameter,num_frames);
    datmat(:,:,1) = repmat(init, 1,N);        
    snaptime_now  = 1;
    
            temp1 = zeros(N,num_particles);
            temp2 = zeros(2,N);
            temp3 = zeros(N,num_particles);
            
    
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

    derivative_alltime = zeros(num_parameters, num_slices); 
 
%%%Initialization%%%
    v = sigV*rnsource;
    initx = repmat(init,1,N);
    data_now = initx;
    timemat = delta*(0:(num_timepts+1));    
    deriv_loglike = zeros(num_parameter,N);
    
%%%Debug : Assessment variables %%%             
    energy = zeros(1, num_slices);  

%%%Monte Carlo Simulation%%%     
    for(m = 2 : (num_timepts+1))
        
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * datmat(1,:,m-1);
        rxn_rate(3,:) = theta(3) * datmat(1,:,m-1);
        rxn_rate(4,:) = theta(4) * datmat(2,:,m-1);
        
        %rxn_rate(1,:) = theta(1); 
        %rxn_rate(2,:) = theta(2) * data_now(1,:);
        %rxn_rate(3,:) = theta(3) * data_now(1,:);
        %rxn_rate(4,:) = theta(4) * data_now(2,:);        
        
        waitbar(m/num_timepts);
        datmathat = datmat(:, :, m-1) + rxn_mat*rxn_rate*delta;
       % datmat(:, :, m)  =  max(datmathat + sigV* sqrt(delta)* rnsource(:,:,m-1),0)  ;
        datmat(:, :, m)  =  datmathat + sigV* sqrt(delta)* rnsource(:,:,m-1);
        
        %%% ftheta part is considered here.
        deriv_loglike(1,:) = deriv_loglike(1,:) +  rxn_mat(:,1)'*(v(:,:,m-1)*sqrt(delta))/(sigV^2);
        deriv_loglike(2,:) = deriv_loglike(2,:) +  data_now(1,:).*(rxn_mat(:,2)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)   ;  
        deriv_loglike(3,:) = deriv_loglike(3,:) +  data_now(1,:).*(rxn_mat(:,3)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)  ;
        deriv_loglike(4,:) = deriv_loglike(4,:) +  data_now(2,:).*(rxn_mat(:,4)'*(v(:,:,m-1)*sqrt(delta)))/(sigV^2)  ;        
        
                
        data_now = data_now + rxn_mat*rxn_rate*delta + sqrt(delta)*v(:,:,m-1);

        
        if(timemat(m) == snaptime(snaptime_now))

            temp1 = repmat(datmat(:,:,m),1, 1, num_particles);
            temp2 = permute(repmat(snapshots(:,:,snaptime_now),1,1,N)...
                ,[1,3,2]);
            temp3 = temp1 - temp2;
            temp3 = min(sum(temp3.^2, 1)/sigW^2, 300);
                        %exp(-735) is the lower limit of precision. 
            tilde_pym =  permute(exp(-temp3), [2,3,1]);             
            sum_tilde_pym = sum(tilde_pym,2);
            sum_tilde_pym = 1./sum_tilde_pym;
            tilde_pym = repmat(sum_tilde_pym,[1,num_particles]) .* tilde_pym;
            deriv_pym = permute(repmat(deriv_loglike, [1,1,num_particles]),[2,3,1]).* ...
                repmat(tilde_pym, [1,1,num_parameter]);
            deriv_py = permute(mean(deriv_pym,1), [3,2,1]) ;   
            tilde_py = repmat(mean(tilde_pym,1),[num_parameter,1]);  
            deriv_pm(:, snaptime_now) = mean(deriv_py./tilde_py,2); 

            tilde_pys(snaptime_now,:) = tilde_py(1,:);
            %scatter(snapshots(1,:,snaptime_now),snapshots(2,:,snaptime_now), 8, tilde_py)

            %%% New Codes %%% 
            xcopy = repmat(data_now,[1,1,num_particles]);             
            ycopy = permute(repmat(snapshots(:,:,snaptime_now),[1,1,N])...
                ,[1,3,2]);
            distance_pair = min(squeeze(sum((xcopy - ycopy).^2,1)/sigW^2), 300);
            p_ymkj = exp(-(distance_pair'));
            p_ymkj = p_ymkj * diag(1./sum(p_ymkj,1));   %tilde_pym in the old version
            tilde_p_ymk = mean(p_ymkj,2);
            dEP_kr = 1/N *p_ymkj * deriv_loglike' ;
            
            derivative_alltime(:,snaptime_now) = compress_snap_wgts(:,snaptime_now)' * ...
                (dEP_kr ./ repmat(tilde_p_ymk,[1, num_parameters]));

            snaptime_now = snaptime_now  +1;

                   if(max(tilde_pym(:)) > min(tilde_pym(:)))
                        include = [include, snaptime_now];
                   end
        end %end of frame
                               

    end %end of time loop 
    derivative = mean(deriv_pm(:,(include-1)),2);
    derivative = mean(derivative_alltime(:,(include-1)),2);
    close(h);
end 

function [rxn_rate] = compute_rxn_rate(theta, data)
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * data(1,:);
        rxn_rate(3,:) = theta(3) * data(1,:);
        rxn_rate(4,:) = theta(4) * data(2,:);
end 
