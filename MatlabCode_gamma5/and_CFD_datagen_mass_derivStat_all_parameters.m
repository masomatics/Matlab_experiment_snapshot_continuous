
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


%         for j = 1:N 
%             deriv_loglike2(1,j) = deriv_loglike2(1,j) + rxn_mat(:,1)'*(datmat(:,j,m) - datmathat(:,j));
%             deriv_loglike2(2,j) = deriv_loglike2(2,j) + datmat(1,j,m-1)* rxn_mat(:,2)'*(datmat(:,j,m) - datmathat(:,j));
%             deriv_loglike2(3,j) = deriv_loglike2(3,j) + datmat(1,j,m-1)* rxn_mat(:,3)'*(datmat(:,j,m) - datmathat(:,j));
%             deriv_loglike2(4,j) = deriv_loglike2(4,j) + datmat(2,j,m-1)* datmat(1,j,m-1)*rxn_mat(:,4)'*(datmat(:,j,m) - datmathat(:,j));
%             deriv_loglike2(5,j) = deriv_loglike2(5,j) + datmat(2,j,m-1)* rxn_mat(:,5)'*(datmat(:,j,m) - datmathat(:,j));
%         end 
        
        
        if(timemat(m) == snaptime(snaptime_now))
        %m
        %snaptime(snaptime_now)
        %    sum_tilde_pym2 = zeros(N, num_particles);

            
% Avoid the for loop, for God's sake.            
%             tic,
%              for(j = 1:N)
%                 sum_tilde_pym2(j) = 0; 
%                 for(k = 1 :num_particles)
%                     tilde_pym2(j,k) = exp(- sum( (snapshots(:,k,snaptime_now) - datmat(:,j,m)).^2) / (sigW^2));
%                     deriv_pym2(j,k,1) = tilde_pym2(j,k) * deriv_loglike(1,j);
%                     deriv_pym2(j,k,2) = tilde_pym2(j,k) * deriv_loglike(2,j);
%                     deriv_pym2(j,k,3) = tilde_pym2(j,k) * deriv_loglike(3,j);
%                     deriv_pym2(j,k,4) = tilde_pym2(j,k) * deriv_loglike(4,j);
%                     deriv_pym2(j,k,5) = tilde_pym2(j,k) * deriv_loglike(5,j);    
%                     
%                     sum_tilde_pym2(j) = sum_tilde_pym2(j) + tilde_pym2(j,k);
%                  end
%                  tilde_pym2(j,:) = tilde_pym2(j,:) / sum_tilde_pym2(j);
%                  deriv_pym2(j,:,:) = deriv_pym2(j,:,:) / sum_tilde_pym2(j);
%              end 
%              deriv_py2 = permute(mean(deriv_pym2,1), [3,2,1]) ;
%              tilde_py2 = mean(tilde_pym2,1);

 
           % tic,
            temp1 = repmat(datmat(:,:,m),1, 1, num_particles);
            temp2 = permute(repmat(snapshots(:,:,snaptime_now),1,1,N)...
                ,[1,3,2]);
            temp3 = temp1 - temp2;
            temp3 = min(sum(temp3.^2, 1)/sigW^2, 300);
                        %exp(-735) is the lower limit of precision. 
            tilde_pym =  permute(exp(-temp3), [2,3,1]);             
%[r,c,v] = ind2sub(size(temp2),find(abs(temp2 - temp1) > 80));
           % toc
 
 
 
 
          %  tic,
            sum_tilde_pym = sum(tilde_pym,2);
            sum_tilde_pym = 1./sum_tilde_pym;
            tilde_pym = repmat(sum_tilde_pym,[1,num_particles]) .* tilde_pym;
            deriv_pym = permute(repmat(deriv_loglike, [1,1,num_particles]),[2,3,1]).* ...
                repmat(tilde_pym, [1,1,num_parameter]);
          %  toc
            

          %  tic,
            deriv_py = permute(mean(deriv_pym,1), [3,2,1]) ;   
            tilde_py = repmat(mean(tilde_pym,1),[num_parameter,1]);  
        
          %  toc
            
            deriv_pm(:, snaptime_now) = mean(deriv_py./tilde_py,2); 
            
            
            
            
            %figure;
            tilde_pys(snaptime_now,:) = tilde_py(1,:);
            %scatter(snapshots(1,:,snaptime_now),snapshots(2,:,snaptime_now), 8, tilde_py)
            
        snaptime_now = snaptime_now + 1 ;  

                    if(max(tilde_pym(:)) > min(tilde_pym(:)))
                        include = [include, snaptime_now];
                    end
        
        end 
                                  
    end
    derivative = mean(deriv_pm(:,(include-1)),2);
    close(h);
end 