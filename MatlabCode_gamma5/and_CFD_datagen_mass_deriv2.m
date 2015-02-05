
%%
% This is the diffusion version of Dave's CFD paper. 
%
% dXt = A Xt dt + dW.  
% 
% rnsource  2 by  sample number by numtimepoints 
% dataset   2 by  sample number by numtimepoints 
% init      2 by 1
%% 

function [timemat, datmat, derivdat, derivdat2] = and_CFD_datagen_mass_deriv2(init, tend, theta, sigV, num_timepts, rnsource, N)
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    
    A = [-theta(3),0; theta(2), -theta(4)]; 
    
    timemat = delta*(0:num_timepts);
    
    datmat = zeros(2, N, num_timepts+1);
    derivdat = NaN(2,N);
    deriv_loglike = zeros(1,N);
    deriv_2nd_loglike = zeros(1,N);
    datmat(:,:,1) = repmat(init, 1,N);    
    matgrowth = repmat([1;0],1,N);
    for(k = 1 : num_timepts)
        
        waitbar(k/num_timepts);
        datmathat = datmat(:, :, k) + theta(1)*delta*matgrowth + (A * datmat(:, :, k))*delta;
        datmat(:, :, k+1)  =  datmathat + sigV* sqrt(delta)* rnsource(:,:,k)  ;
        
        
        deriv_loglike = deriv_loglike +   sum( (datmat(:, :, k+1)- datmathat).*([-1, 0 ;0,0] * datmat(:, :, k)), 1) /(sigV^2) ;
        deriv_2nd_loglike = deriv_2nd_loglike + (-1)*  sum(([-1, 0 ;0,0] * datmat(:, :, k)).^2, 1) /(sigV^2)*delta; 
    end 
    derivdat =datmat(:,:,num_timepts+1)*diag(deriv_loglike);
   % derivdat2= datmat(:,:,num_timepts+1)*(diag(deriv_loglike.^2));
    derivdat2= datmat(:,:,num_timepts+1)*(diag(deriv_loglike.^2 +deriv_2nd_loglike ))  ;
    close(h);
end 