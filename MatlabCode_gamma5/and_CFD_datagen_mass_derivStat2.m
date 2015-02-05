
%%
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
%% 

function [datmat, tilde_pys, derivative] = and_CFD_datagen_mass_derivStat2...
                        (init, tend, theta, sigV, sigW, num_timepts, rnsource, snapshots, snaptime, N)
    
    [species, num_particles, num_frames] = size(snapshots);  
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    
    A = [-theta(1),0; theta(2), -1]; 
    
    timemat = delta*(0:num_timepts);
    
    datmat = zeros(2, N, num_timepts +1);
    
    %\beta_m^k
    derivdat = NaN(num_particles,num_frames);
    
    deriv_loglike1 = zeros(1,N);
    deriv_loglike2 = zeros(1,N);

    deriv_pym1 = zeros(N, num_particles);
    deriv_pym2 = zeros(N, num_particles);

    tilde_pym = zeros(N, num_particles);
    sum_tilde_pym = zeros(N,1);
    
    deriv_py1 = zeros(1,num_particles);
    deriv_py2 = zeros(1,num_particles);

    tilde_py = zeros(1,num_particles);
    tilde_pys= zeros(num_frames, num_particles);
    
    deriv_pm = zeros(1,num_frames);
    
    
    datmat(:,:,1) = repmat(init, 1,N);    
    matgrowth = repmat([1;0],1,N);
    
    snaptime_now  = 1;
    
    
            temp1 = zeros(N,num_particles);
            temp2 = zeros(2,N);
            temp3 = zeros(N,num_particles);
            
    
    for(m = 2 : num_timepts+1)
        
        waitbar(m/num_timepts);
        datmathat = datmat(:, :, m-1) + 2*delta*matgrowth + (A * datmat(:, :, m-1))*delta;
        datmat(:, :, m)  =  datmathat + sigV* sqrt(delta)* rnsource(:,:,m-1)  ;
        
        
        deriv_loglike1 = deriv_loglike1 + sum( (datmat(:, :, m)- ...
             datmathat).*([-1, 0 ;0,0] * datmat(:, :, m-1)), 1) /(sigV^2) ;
        deriv_loglike2 = deriv_loglike2 + sum( (datmat(:, :, m)- ...
             datmathat).*([0, 0 ;1,0] * datmat(:, :, m-1)), 1) /(sigV^2) ;        
        if(timemat(m) == snaptime(snaptime_now))
        %m
            %sum_tilde_pym = zeros(N, num_particles);

            
% Avoid the for loop, for God's sake.            
%             tic,
%             for(k = 1 :num_particles)
% 
%                 for(j = 1:N)
%                    %% why is this part so slow 
%                    tilde_pym(j,k) = exp(- sum( (snapshots(:,k,snaptime_now) - datmat(:,j,m)).^2) / (sigW^2));
%                    deriv_pym(j,k) = tilde_pym(j,k) * deriv_loglike(j);
%                 %   sum_tilde_pym(j) = sum_tilde_pym(j) + tilde_pym(k,j);
%                 end
%             end 
%             toc
            

%             tic, 
%             for(k = 1 : num_particles)
%                 temp1 = repmat(snapshots(:,k,snaptime_now), 1, N);
%                 temp2 = temp1 - datmat(:,:,m);
%                 temp3(:,k) = sum(temp2.^2, 1)/sigW^2;
%             end
%             tilde_pym = exp(-temp3);
%             deriv_pym = diag(deriv_loglike)*tilde_pym ;
%             toc
 %           display(tilde_pym2);

            
 
           % tic,
            temp1 = repmat(datmat(:,:,m),1, 1, num_particles);
            temp2 = permute(repmat(snapshots(:,:,snaptime_now),1,1,N)...
                ,[1,3,2]);
            temp3 = temp1 - temp2;
            temp3 = min(sum(temp3.^2, 1)/sigW^2, 500);
                        %exp(-735) is the lower limit of precision. 
            tilde_pym =  permute(exp(-temp3), [2,3,1]);             
          %  deriv_pym1 = diag(deriv_loglike1)*tilde_pym ;
          %  deriv_pym2 = diag(deriv_loglike2)*tilde_pym ;

           % toc
 
 
 
 
          %  tic,
            sum_tilde_pym = sum(tilde_pym,2);
            sum_tilde_pym = 1./sum_tilde_pym;
          %  tic,
          %  deriv_pym1 = diag(deriv_loglike1) * tilde_pym ;
          %  deriv_pym2 = diag(deriv_loglike2) * tilde_pym ;
          %  toc
          %  tic,
            deriv_pym1 = repmat(deriv_loglike1', 1,num_particles).* tilde_pym;
            deriv_pym2 = repmat(deriv_loglike2', 1,num_particles).* tilde_pym;
          %  toc
          %tilde_pym = diag(sum_tilde_pym)* tilde_pym;
            tilde_pym = repmat(sum_tilde_pym,1,num_particles).*tilde_pym;

          %  deriv_pym1 = diag(sum_tilde_pym)* deriv_pym1;
          %  deriv_pym2 = diag(sum_tilde_pym)* deriv_pym2;

          deriv_pym1= repmat(sum_tilde_pym,1,num_particles).*deriv_pym1;
          deriv_pym2= repmat(sum_tilde_pym,1,num_particles).*deriv_pym2;
          %  toc
            

          %  tic,
            deriv_py1 = mean(deriv_pym1,1);   
            deriv_py2 = mean(deriv_pym2,1);   

            tilde_py = mean(tilde_pym,1);  
        
          %  toc
            
            deriv_pm1(snaptime_now) = mean(deriv_py1./tilde_py); 
            deriv_pm2(snaptime_now) = mean(deriv_py2./tilde_py); 
            
            %figure;
            tilde_pys(snaptime_now,:) = tilde_py;
            %scatter(snapshots(1,:,snaptime_now),snapshots(2,:,snaptime_now), 8, tilde_py)
            
        snaptime_now = snaptime_now + 1 ;  

        end 
                                  
    end
    
    derivative = [mean(deriv_pm1),mean(deriv_pm2)] ;
    close(h);
end 