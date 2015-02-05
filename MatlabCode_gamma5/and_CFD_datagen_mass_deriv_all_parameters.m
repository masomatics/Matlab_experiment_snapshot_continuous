
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

function [timemat, datmat, derivdat] = and_CFD_datagen_mass_derivStat_test...
                       (init, tend, theta, sigV, num_timepts, rnsource, N)
    
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    num_parameter = 4;
    
    rxn_mat =   [1 0 -1 0; 0 1 0 -1]; 
    rxn_rate = zeros(4,N);
    
    timemat = delta*(0:(num_timepts+1));
    
    datmat = zeros(2, N, num_timepts +1);
    
    %\beta_m^k
    derivdat = NaN(2, N, num_parameter);
    
    deriv_loglike = zeros(num_parameter,N);
   % deriv_loglike2 = zeros(5,N);
    
    datmat(:,:,1) = repmat(init, 1,N);    
       
    
            
    for(m = 2 : (num_timepts+1))
        
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * datmat(1,:,m-1);
        rxn_rate(3,:) = theta(3) * datmat(1,:,m-1);
        rxn_rate(4,:) = theta(4) * datmat(2,:,m-1);
        
        
        
        waitbar(m/num_timepts);
        datmathat = datmat(:, :, m-1) + rxn_mat*rxn_rate*delta;
       % datmat(:, :, m)  =  max(datmathat + sigV* sqrt(delta)* rnsource(:,:,m-1),0)  ;
        datmat(:, :, m)  =  datmathat + sigV* sqrt(delta)* rnsource(:,:,m-1);
        
        
        %deriv_loglike(1,:) = deriv_loglike(1,:) +   sum( (datmat(:, :, k+1)- datmathat).*([-1, 0 ;0,0] * datmat(:, :, k)), 1) /(sigV^2) ;
        
        deriv_loglike(1,:) = deriv_loglike(1,:) +  rxn_mat(:,1)'*(datmat(:,:,m) - datmathat)/(sigV^2);
        deriv_loglike(2,:) = deriv_loglike(2,:) +  datmat(1,:,m-1).*(rxn_mat(:,2)'*(datmat(:,:,m) - datmathat))/(sigV^2)   ;  
        deriv_loglike(3,:) = deriv_loglike(3,:) +  datmat(1,:,m-1).*(rxn_mat(:,3)'*(datmat(:,:,m) - datmathat))/(sigV^2)  ;
        deriv_loglike(4,:) = deriv_loglike(4,:) +  datmat(2,:,m-1).*(rxn_mat(:,4)'*(datmat(:,:,m) - datmathat))/(sigV^2)  ;

    end
    derivdat(:,:,1) = datmat(:,:,num_timepts+1).* repmat(deriv_loglike(1,:),[2,1]);
    derivdat(:,:,2) = datmat(:,:,num_timepts+1).* repmat(deriv_loglike(2,:),[2,1]);
    derivdat(:,:,3) = datmat(:,:,num_timepts+1).* repmat(deriv_loglike(3,:),[2,1]);
    derivdat(:,:,4) = datmat(:,:,num_timepts+1).* repmat(deriv_loglike(4,:),[2,1]);
    
    
    close(h);
end 