
%
% Solution to ther Anderson CFD. Use this to confirm the formulas for the
% first and the second derivative. 
%
%
theta(1) = 2 ; 
theta(2) = 10;
theta(3) = 1/4;
theta(4) = 1;
t = 30; 
X0 = [0; 0]

A = [-theta(3), 0 ; theta(2) , -theta(4)]; 
b = [theta(1); 0];

P   = [theta(4) - theta(3) , 0; theta(2), 1];
Pin = [1, 0 ; -theta(2), theta(4) - theta(3)]./(theta(4) - theta(3)) ;




D = diag([-theta(3) , - theta(4)]);
expDt = diag([ exp(-theta(3)*t) , exp(- theta(4)*t)])
intD_t = diag([(exp(theta(3)*t)-1)/theta(3), (exp(theta(4)*t)-1)/theta(4)]) 

expAt = P * expDt * Pin;

int_expA_t_s =   P*(expDt * intD_t)*Pin

Xt = expAt * X0 + int_expA_t_s * b

%Full Formula
Xt = P* diag([exp(-theta(3)*t), exp(-theta(4)*t)]) * Pin  * X0+ ... 
P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*Pin...
*[theta(1);0]


DD = diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))]);
rho_3P = [-1, 0; 0,0];
rho_3Pin = [0, 0 ; 0, -1] ./(theta(4) - theta(3)) + [1, 0 ; -theta(2), theta(4) - theta(3)]./((theta(4) - theta(3))^2)
rho_3DD = diag([-1/(theta(3)^2)*(1- exp(-theta(3)*t)) + 1/theta(3) * t* exp(-theta(3)*t),0]);


rho_3Xt = (P* diag([-t* exp(-theta(3)*t), 0 ]) * Pin + ... 
rho_3P*     diag([exp(-theta(3)*t), 0 ]) *Pin + ...
P*     diag([exp(-theta(3)*t), 0 ]) *rho_3Pin) * X0+ ... 
(P*diag([-1/(theta(3)^2)*(1- exp(-theta(3)*t)) + 1/theta(3) * t* exp(-theta(3)*t),0])*Pin +...
rho_3P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*Pin + ...
P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*rho_3Pin)...
*[theta(1);0]

rho_3Xt = (P*rho_3DD*Pin + rho_3P*DD*Pin + P * DD * rho_3Pin)* [theta(1);0]


%I will ignore the first term since X0 in our example is the zero vector. 
rho_33P = [0, 0; 0,0];
rho_33Pin = [0, 0 ; 0, -1] ./(theta(4) - theta(3))^2 + ...
    [0, 0 ; 0, -1]./((theta(4) - theta(3))^2) + ...
    [1, 0 ; -theta(2), theta(4) - theta(3)]./((theta(4) - theta(3))^3)*2 ;   

rho_33DD = diag([ ...
2/(theta(3)^3)*(1- exp(-theta(3)*t))-...    
2/(theta(3)^2)*t*exp(- theta(3)*t)+...
1/(theta(3))* (-t^2)*exp(- theta(3)*t) ...
,0]);

rho_33Xt= ( rho_3P*rho_3DD*Pin + P*rho_33DD*Pin + P*rho_3DD*rho_3Pin + ...
  rho_33P*DD*Pin + rho_3P*rho_3DD*Pin  + rho_3P*DD*rho_3Pin + ...
  rho_3P*DD*rho_3Pin + P*rho_3DD*rho_3Pin + P*DD*rho_33Pin)* [theta(1);0] 

%%
initX = X0;
sigV = 2;
sigW = 2;
%num_timepts = 2500;
num_timepts = 500;
%Ntry = 10000;
Ntry = 200
tic,
rnsource = randn([2, Ntry, num_timepts]);
toc
tic,

%% two parameters version
tendtest = t
initx
theta
[sigV,sigW]

[timepts,datapts, derivpts, derivpts2] = and_CFD_datagen_mass_deriv2(initx, tendtest, theta, sigV, num_timepts, rnsource, Ntry);

display(['Our estimate of mean1 is: ', num2str(mean(datapts(1,:,num_timepts))), ' real value is ', num2str(Xt(1))  ]);  
display(['Our estimate of mean2 is: ', num2str(mean(datapts(2,:,num_timepts))), ' real value is ', num2str(Xt(2)) ]); 

display(['Our estimate of derivmean1 is: ', num2str(mean(derivpts(1,:))), ' real value is ', num2str(rho_3Xt(1)) ]); 
display(['the variance of the esimtate of derivmean2 is: ',...
    num2str(var(derivpts(1,:))/Ntry)]); 

display(['Our estimate of derivmean2 is: ', num2str(mean(derivpts(2,:))), ' real value is ', num2str(rho_3Xt(2))]); 
display(['the variance of the esimtate of derivmean2 is: ',...
    num2str(var(derivpts(2,:))/Ntry)]); 

display(['Our estimate of 2nd derivmean1 is: ', num2str(mean(derivpts2(1,:))), ' real value is ', num2str(rho_33Xt(1))]); 
display(['the variance of the esimtate of 2nd derivmean2 is: ',...
    num2str(var(derivpts2(1,:))/Ntry)]); 

display(['Our estimate of 2nd derivmean2 is: ', num2str(mean(derivpts2(2,:))), ' real value is ', num2str(rho_33Xt(2))]); 
display(['the variance of the esimtate of 2nd derivmean2 is: ',...
    num2str(var(derivpts2(2,:))/Ntry)]); 

toc

%% ALL parameters version ... upto first derivatives. 


[timepts,datapts_all_p, derivpts_all_p] = and_CFD_datagen_mass_deriv_all_parameters(initx, tendtest, theta, sigV, num_timepts, rnsource, Ntry);

display(['Our estimate of mean1 is: ', num2str(mean(datapts_all_p(1,:,num_timepts))), ' real value is ', num2str(Xt(1))  ]);  
display(['Our estimate of mean2 is: ', num2str(mean(datapts_all_p(2,:,num_timepts))), ' real value is ', num2str(Xt(2)) ]); 

display(['Our estimate of derivmean1 is: ', num2str(mean(derivpts_all_p(1,:,3))), ' real value is ', num2str(rho_3Xt(1)) ]); 
display(['the variance of the esimtate of derivmean2 is: ',...
    num2str(var(derivpts_all_p(1,:))/Ntry)]); 

display(['Our estimate of derivmean2 is: ', num2str(mean(derivpts_all_p(2,:,3))), ' real value is ', num2str(rho_3Xt(2))]); 
display(['the variance of the esimtate of derivmean2 is: ',...
    num2str(var(derivpts_all_p(2,:))/Ntry)]); 


%% All parameters version .. upto second derivatives. This part is not done yet.
display(['Our estimate of 2nd derivmean1 is: ', num2str(mean(derivpts2(1,:))), ' real value is ', num2str(rho_33Xt(1))]); 
display(['the variance of the esimtate of 2nd derivmean2 is: ',...
    num2str(var(derivpts2(1,:))/Ntry)]); 

display(['Our estimate of 2nd derivmean2 is: ', num2str(mean(derivpts2(2,:))), ' real value is ', num2str(rho_33Xt(2))]); 
display(['the variance of the esimtate of 2nd derivmean2 is: ',...
    num2str(var(derivpts2(2,:))/Ntry)]); 

toc

% [sigV,sigW]
% 
% ans =
% 
%      2     2
% 
% initx =
% 
%      1
%      0
% 
% theta =
% 
%     2.0000   10.0000    0.3333    1.0000
% 
% 
% Ntry =
% 
%        80000
% 
% Elapsed time is 6.538397 seconds.
% 
% tendtest =
% 
%     30
% 
% Our estimate of mean1 is: 5.9885 real value is 5.9997
% Our estimate of mean2 is: 59.8894 real value is 59.9959
% Our estimate of derivmean1 is: -17.7275 real value is -17.991
% the variance of the esimtate of derivmean2 is: 0.17618
% Our estimate of derivmean2 is: -176.725 real value is -179.8713
% the variance of the esimtate of derivmean2 is: 17.2864
% Our estimate of 2nd derivmean1 is: 114.4608 real value is 107.7009
% the variance of the esimtate of 2nd derivmean2 is: 195.0409
% Our estimate of 2nd derivmean2 is: 1132.1977 real value is 1075.8997
% the variance of the esimtate of 2nd derivmean2 is: 19280.963
% Elapsed time is 265.769487 seconds.
% 
