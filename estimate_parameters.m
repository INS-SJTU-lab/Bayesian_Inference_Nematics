
clear;clc;close all

x_scale = 1329.9/1024;  % 1 pixel = 1.3 um
t_scale = 360; % 1 frame = 360s;


center= [386,386];
%the size of picture(pixel)
LX = 771;
LY = 771;
dis = 32;
xlist1 = center(1):dis:LX;
xlist2 = fliplr(center(1):-dis:1);
xlist = [xlist2(1:end-1),xlist1];
ylist1 = center(2):dis:LY;
ylist2 = fliplr(center(2):-dis:1);
ylist = [ylist2(1:end-1),ylist1];
[xmesh,ymesh] = meshgrid(xlist,ylist);

LXv = 97;
LYv = 97;
center= [floor(LXv/2)+1,floor(LXv/2)+1];
%the size of picture(pixel)
dis = 4;
xlist1 = center(1):dis:LXv;
xlist2 = fliplr(center(1):-dis:1);
xlistv = [xlist2(1:end-1),xlist1];
ylist1 = center(2):dis:LYv;
ylist2 = fliplr(center(2):-dis:1);
ylistv = [ylist2(1:end-1),ylist1];
[xmeshv,ymeshv] = meshgrid(xlistv,ylistv);


for kk = 1:1:4 % different regions

    % the average of velocity field
    load('velocity_field.mat');
    VX_ave = sum(V_grid_t{kk,1},3)./sum(V_grid_t{kk,3},3)/60; 
    VY_ave = sum(V_grid_t{kk,2},3)./sum(V_grid_t{kk,3},3)/60;
    sigma1 = 4;
    VX_ave = imgaussfilt(VX_ave, sigma1);  % filtered
    VY_ave = imgaussfilt(VY_ave, sigma1);  % filtered
    sVX_ave = VX_ave(ylistv,xlistv);  % Sampling every 2 lattice points
    sVY_ave = VY_ave(ylistv,xlistv);  % Sampling every 2 lattice points
    sV_amplitude = sqrt(sVX_ave.^2+sVY_ave.^2);   

    load('densityfield.mat');
    sigma = 32;
    rho = densityfield_ave(:,:,kk);
    rho = imgaussfilt(rho, sigma);
    % rho = rho/normal;   % Normalize the density by the maximum value
    normal = mean(rho(:));
    rho = rho/normal;   % Normalize the density by the maximum value 
    [dxrho,dyrho] = gradient(rho,x_scale);  % partial_x QPxx
    sdxrho = dxrho(ylist,xlist);
    sdyrho = dyrho(ylist,xlist);
    srho = rho(ylist,xlist);

    % the average of Q
    load('orientation_field.mat');
    Qxx = orientation_grid_t{kk,1};
    Qxy = orientation_grid_t{kk,2};
    S_ave = sqrt(Qxx.^2 + Qxy.^2);   
    sQxx = Qxx(ylist,xlist);
    sQxy = Qxy(ylist,xlist); 

    rhoQxx = rho.*Qxx;
    [dxrhoQxx,dyrhoQxx] = gradient(rhoQxx,x_scale);
    rhoQxy = rho.*Qxy;
    [dxrhoQxy,dyrhoQxy] = gradient(rhoQxy,x_scale);
    dxrhoQxx = imgaussfilt(dxrhoQxx, sigma);
    dyrhoQxx = imgaussfilt(dyrhoQxx, sigma);
    sdxrhoQxx = dxrhoQxx(ylist,xlist);
    sdyrhoQxx = dyrhoQxx(ylist,xlist);
    dxrhoQxy = imgaussfilt(dxrhoQxy, sigma);
    dyrhoQxy = imgaussfilt(dyrhoQxy, sigma);
    sdxrhoQxy = dxrhoQxy(ylist,xlist);
    sdyrhoQxy = dyrhoQxy(ylist,xlist);
    srhoQxx = rhoQxx(ylist,xlist);
    srhoQxy = rhoQxy(ylist,xlist);

    % force field
    % pi2 term   
    Fx1 = sdxrhoQxx + sdyrhoQxy;
    Fy1 = sdxrhoQxy - sdyrhoQxx;

    % gamma2 term
    Fx2 = srho.*sQxx.*(sdxrhoQxx+sdyrhoQxy)+srho.*sQxy.*(sdxrhoQxy-sdyrhoQxx);
    Fy2 = srho.*sQxy.*(sdxrhoQxx+sdyrhoQxy)+srho.*sQxx.*(-sdxrhoQxy+sdyrhoQxx);

    % gamma1 term
    Fx3 = srho.*sQxx.*(sdxrhoQxx-sdyrhoQxy)+srho.*sQxy.*(sdxrhoQxy+sdyrhoQxx);
    Fy3 = -srho.*sQxy.*(sdxrhoQxx-sdyrhoQxy)+ srho.*sQxx.*(sdxrhoQxy+sdyrhoQxx);

    % friction term
    frix = srho.*sQxx.*srho.*sVX_ave + srho.*sQxy.*srho.*sVY_ave;
    friy = srho.*sQxy.*srho.*sVX_ave - srho.*sQxx.*srho.*sVY_ave;

    % chiral term
    chiralx = -srho.*sVY_ave;
    chiraly = srho.*sVX_ave;

    % diffusion term
    diffx = -sdxrho;
    diffy = -sdyrho;

    % anidiffusion term    
    anidiffx = -(srho.*sQxx.*sdxrho+srho.*sQxy.*sdyrho);
    anidiffy = -(srho.*sQxy.*sdxrho-srho.*sQxx.*sdyrho);

    % common meshgrid, exclude region near the core and near the edge
    center= [386,386];
    %the size of picture(pixel)
    LX = 771;
    LY = 771;
    dis = 32;
    xlist1 = center(1):dis:LX;
    xlist2 = fliplr(center(1):-dis:1);
    xlist = [xlist2(1:end-1),xlist1];
    ylist1 = center(2):dis:LY;
    ylist2 = fliplr(center(2):-dis:1);
    ylist = [ylist2(1:end-1),ylist1];
    [xmesh,ymesh] = meshgrid(xlist,ylist);
    DX = xmesh - center(1);
    DY = ymesh - center(2);
    ddr_mat = sqrt(DX.^2+DY.^2);
    radius2 = 450; %看看这选多少合适
    radius2 = radius2*1024/1329.9; %unit in pixel
    radius1 = 100;
    radius1 = radius1*1024/1329.9; %unit in pixel
    validid = ddr_mat>=radius1 & ddr_mat<=radius2;
    phi_mat = atan2(DY,DX);
    
    % average angle
    theta = atan2(sQxy, sQxx)/2;   % Angle from Q11 and Q12   
    % dot = cos(theta).*cos(phi_mat) + sin(theta).*sin(phi_mat);
    % theta(dot<0) = theta(dot<0) + pi;
    [F1para_s,F1perp_s] = project2ori(Fx1,Fy1,theta); 
    [F2para_s,F2perp_s] = project2ori(Fx2,Fy2,theta); 
    [F3para_s,F3perp_s] = project2ori(Fx3,Fy3,theta); 
    [fripara_s,friperp_s] = project2ori(frix,friy,theta); 
    [chiralpara_s,chiralperp_s] = project2ori(chiralx,chiraly,theta);
    [diffpara_s,diffperp_s] = project2ori(diffx,diffy,theta);
    [anidiffpara_s,anidiffperp_s] = project2ori(anidiffx,anidiffy,theta);
    % project V to Q direction
    [Vpara_s,Vperp_s] = project2ori(sVX_ave,sVY_ave,theta);


    Fx_total = [Fx1(validid),Fx2(validid),Fx3(validid),diffx(validid),anidiffx(validid)];
    Vx_total = [srho(validid).*sVX_ave(validid),frix(validid),chiralx(validid)];  
    Fy_total = [Fy1(validid),Fy2(validid),Fy3(validid),diffy(validid),anidiffy(validid)];
    Vy_total = [srho(validid).*sVY_ave(validid),friy(validid),chiraly(validid)];  
    Fx_total_region(:,:,kk) = Fx_total;
    Vx_total_region(:,:,kk) = Vx_total;
    Fy_total_region(:,:,kk) = Fy_total;
    Vy_total_region(:,:,kk) = Vy_total;
    theta_total_region(:,:,kk) = theta(validid);
    sQxx_total_region(:,:,kk) = srhoQxx(validid);
    sQxy_total_region(:,:,kk) = srhoQxy(validid);
      
    F_total = [[F1para_s(validid);F1perp_s(validid)],[F2para_s(validid);F2perp_s(validid)],...
        [F3para_s(validid);F3perp_s(validid)],[diffpara_s(validid);diffperp_s(validid)],[anidiffpara_s(validid);anidiffperp_s(validid)]];
    V_total = [srho(validid).*Vpara_s(validid);srho(validid).*Vperp_s(validid)];
    Fpara_total_region(:,:,kk) = [F1para_s(validid),F2para_s(validid),F3para_s(validid),diffpara_s(validid),anidiffpara_s(validid)];
    Vpara_total_region(:,:,kk) = srho(validid).*Vpara_s(validid);
    Fperp_total_region(:,:,kk) = [F1perp_s(validid),F2perp_s(validid),F3perp_s(validid),diffperp_s(validid),anidiffperp_s(validid)];
    Vperp_total_region(:,:,kk) = srho(validid).*Vperp_s(validid);
    F_total_region_pp(:,:,kk) = F_total;
    V_total_region_pp(:,:,kk) = V_total;    
end
x_total = xmesh(validid);y_total = ymesh(validid);
save('H06_velocity_sano_2rho_second.mat','Fpara_total_region','Vpara_total_region','Fperp_total_region','Vperp_total_region',...
    'Fx_total_region','Vx_total_region','Fy_total_region','Vy_total_region','theta_total_region','F_total_region_pp','V_total_region_pp',...
    'x_total','y_total','sQxx_total_region','sQxy_total_region');


%% 
epsilon_total = zeros(1,size(Fx_total_region,3));
R2_total = zeros(1,size(Fx_total_region,3));
coeff_total = zeros(size(Fx_total_region,2),size(Fx_total_region,3));

for kk = 1:size(Fx_total_region,3)

    epsilon = 0:0.1:1;
    R2_matrix = zeros(1,length(epsilon));
    coeff_matrix = zeros(size(Fx_total_region,2),length(epsilon));
    error_matrix = R2_matrix;
    for ii = 1:length(epsilon)
        Vx = Vx_total_region(:,1,kk) - epsilon(ii)*Vx_total_region(:,2,kk);
        Vy = Vy_total_region(:,1,kk) - epsilon(ii)*Vy_total_region(:,2,kk);
        Fx = Fx_total_region(:,:,kk);
        Fy = Fy_total_region(:,:,kk);
        F_total = [Fx;Fy];
        V_total = [Vx;Vy];
        % coeff = lsqr(F_total,V_total,1e-6,20);
        [coeff,error] = estimatepara(F_total,V_total);   

        Qxx = sQxx_total_region(:,:,kk);
        Qxy = sQxy_total_region(:,:,kk);
        Vx_predict = (Fx*coeff.*(1+epsilon(ii)*Qxx)+Fy*coeff.*(epsilon(ii)*Qxy))./(1-epsilon(ii)^2*(Qxy.^2+Qxx.^2));
        Vy_predict = (Fy*coeff.*(1-epsilon(ii)*Qxx)+Fx*coeff.*(epsilon(ii)*Qxy))./(1-epsilon(ii)^2*(Qxy.^2+Qxx.^2));    
        V_predict = [Vx_predict;Vy_predict];
        Vx_actual = Vx_total_region(:,1,kk);
        Vy_actual = Vy_total_region(:,1,kk);
        V_actual = [Vx_actual;Vy_actual];
        R2 = 1-sum((V_actual-V_predict).^2)/sum((V_actual-mean(V_actual)).^2);
        R2_matrix(ii) = R2;  
        coeff_matrix(:,ii)=coeff;
        error_matrix(:,ii) = error;
    end
    sort_R2_matrix = sort(R2_matrix);             
    epsilon_index = find(R2_matrix==sort_R2_matrix(end));
    epsilon_max = epsilon(epsilon_index);
    if sort_R2_matrix(end)-sort_R2_matrix(end-1)<=0.00001
        break;
    end

    for gap = [0.01,0.001,0.0001]
        epsilon=max(0,(epsilon_max-gap*10)):gap:min((epsilon_max+gap*10),1);
        R2_matrix = zeros(1,length(epsilon));
        coeff_matrix = zeros(size(Fx_total_region,2),length(epsilon));
        error_matrix = R2_matrix;
        for ii = 1:length(epsilon)
            Vx = Vx_total_region(:,1,kk) - epsilon(ii)*Vx_total_region(:,2,kk);
            Vy = Vy_total_region(:,1,kk) - epsilon(ii)*Vy_total_region(:,2,kk);
            Fx = Fx_total_region(:,:,kk);
            Fy = Fy_total_region(:,:,kk);
            F_total = [Fx;Fy];
            V_total = [Vx;Vy];
            % coeff = lsqr(F_total,V_total,1e-6,20);
            [coeff,error] = estimatepara(F_total,V_total);

            Qxx = sQxx_total_region(:,:,kk);
            Qxy = sQxy_total_region(:,:,kk);
            Vx_predict = (Fx*coeff.*(1+epsilon(ii)*Qxx)+Fy*coeff.*(epsilon(ii)*Qxy))./(1-epsilon(ii)^2*(Qxy.^2+Qxx.^2));
            Vy_predict = (Fy*coeff.*(1-epsilon(ii)*Qxx)+Fx*coeff.*(epsilon(ii)*Qxy))./(1-epsilon(ii)^2*(Qxy.^2+Qxx.^2));
            V_predict = [Vx_predict;Vy_predict];
            Vx_actual = Vx_total_region(:,1,kk);
            Vy_actual = Vy_total_region(:,1,kk);
            V_actual = [Vx_actual;Vy_actual];
            R2 = 1-sum((V_actual-V_predict).^2)/sum((V_actual-mean(V_actual)).^2);
            R2_matrix(ii) = R2;
            coeff_matrix(:,ii) = coeff;
            error_matrix(:,ii) = error;
        end
        sort_R2_matrix = sort(R2_matrix);
        epsilon_index = find(R2_matrix==sort_R2_matrix(end));
        epsilon_max = epsilon(epsilon_index);
        if sort_R2_matrix(end)-sort_R2_matrix(end-1)<=0.00001
            break;
        end
    end
    epsilon = epsilon_max;
    coeff = coeff_matrix(:,epsilon_index);
    R2 = sort_R2_matrix(end);

    epsilon_total(kk) = epsilon;
    R2_total(kk) = R2;
    coeff_total(:,kk) = coeff;
end


function [Fpara,Fperp] = project2ori(Fx,Fy,theta)
cosQ = cos(theta);
sinQ = sin(theta);

Fpara =  Fx.*cosQ + Fy.*sinQ;  % dot(nablaQ, Q_para) 
Fperp = -Fx.*sinQ + Fy.*cosQ;  % dot(nablaQ, Q_perp)
end

function [coeff_Bayesian,error,ridgepara,coeff_Bayesian_changing] = estimatepara(F_total,V_total)

s2 = 0.000001; s02=0.0001;
s_pi0 = 0.0001; % estimate of the variance of pi
pi0 = 0.16;

coeff_Bayesian_changing = [];
Lambda_changing = [];
for i = 1:40

    Lambda = s2/s02;   
    Lambda_pi0 = s2/s_pi0;
    XX = F_total.'*F_total;
    XY = F_total.'*V_total;
    XX(4,4)=XX(4,4)+(Lambda_pi0-Lambda);
    XY(4)=XY(4)+Lambda_pi0*pi0;
    B2 = inv(XX + Lambda*eye(4))*XY;
    Y = V_total;
    X = F_total;
    Predict = X*B2; % Prediction of (uf, vf)
    normer = (Y - Predict)'*(Y - Predict);
    s2 = normer/length(Y);
    s02 = B2([1:3])'*B2([1:3])/3;
    Lambda_n = s2/s02;
    coeff_Bayesian_changing=[coeff_Bayesian_changing,B2];
    Lambda_changing = [Lambda_changing,Lambda];
    if abs(Lambda - Lambda_n)/Lambda <= 0.001
        break
    end
end
coeff_Bayesian = B2;
error = abs(Lambda - Lambda_n)/Lambda;
ridgepara = Lambda;
end