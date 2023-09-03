function [Mz_Fycell] = Mz_lateralforce_combined(Xbestcell,alpha,F_z,p_i,gamma,kappa)
%% Inputs

gamma = gamma*0.0174533; %degrees to radians
alpha = alpha*0.0174533; %degrees to radians

F_z0 = round(F_z/25)*25; %nominal load
p_i0 = round(p_i); %nominal pressure

df_z = (F_z- F_z0)./F_z;
dp_i = (p_i-p_i0)./p_i;

%% Parameters

[p_cy1,p_dy1,p_dy2,p_dy3,p_ey1,p_ey2,p_ey3,p_ey4,p_ey5,p_hy1,p_hy2,p_ky1,p_ky2,p_ky3,...
    p_ky4,p_ky5,p_ky6,p_ky7,p_py1,p_py2,p_py3,p_py4,p_py5,p_vy1,p_vy2,p_vy3,p_vy4,r_by1,r_by2,...
    r_by3,r_by4,r_cy1,r_ey1,r_ey2,r_hy1,r_hy2,r_vy1,r_vy2,r_vy3,r_vy4,r_vy5,r_vy6] = Xbestcell{:};

lambda_cy = 1;       % shape factor   
lambda_ey = 1;       % curvature 
lambda_hy = 1;       % horizontal shift
lambda_kyalpha = 1;  % cornering stiffness
lambda_kygamma = 1;  % camber force stiffness
lambda_muy = 1;      % peak friction coefficient
lambda_vy = 1;       % vertical shift

lambda_ykappa = 1;   % influence on F_y(alpha)
lambda_vykappa = 1;  % induced ply-steer F_y

%%  Magic Formula Equations
K_yalpha = p_ky1.*F_z0.*(1+p_py1.*dp_i).*sin(p_ky4.*atan(F_z./...
    ((p_ky2 + p_ky5*gamma.^2).*(1+p_py2*dp_i).*F_z0))).*...
    (1-p_ky3.*abs(gamma)).*lambda_kyalpha;
K_ygamma = (p_ky6+p_ky7.*df_z).*(1+p_py5.*dp_i).*F_z.*lambda_kygamma;
S_vy0 = F_z.*(p_vy1+p_vy2.*df_z).*lambda_vy.*lambda_muy;
S_vygamma = F_z.*(p_vy3+p_vy4.*df_z).*gamma.*lambda_kygamma.*lambda_muy;
S_vy = S_vy0+S_vygamma;
S_hy0 = (p_hy1+p_hy2.*df_z).*lambda_hy;
S_hygamma = (K_ygamma.*gamma-S_vygamma)./K_yalpha;
S_hy = S_hy0+S_hygamma;
alpha_y = alpha+S_hy;
mu_y = (p_dy1+p_dy2.*df_z).*(1-p_dy3.*gamma.^2).*...
    (1+p_py3.*dp_i+p_py4.*dp_i.^2).*lambda_muy;
C_y = p_cy1.*lambda_cy;
D_y = mu_y.*F_z;
E_y = (p_ey1+p_ey2.*df_z).*(1+p_ey5.*gamma.^2-(p_ey3+p_ey4.*gamma).*sign(alpha_y)).*lambda_ey;
B_y = K_yalpha./C_y./D_y;

%% Combined Slip
B_ykappa = (r_by1 + r_by4.*gamma.^2).*cos(atan(r_by2.*(alpha - r_by3))).*lambda_ykappa;
C_ykappa = r_cy1;
E_ykappa = r_ey1 + r_ey2.*df_z;
S_hykappa = r_hy1 + r_hy2.*df_z; 
kappa_s = kappa + S_hykappa;
D_vykappa = mu_y.*F_z.*(r_vy1 + r_vy2.*df_z + r_vy3.*gamma).*cos(atan(r_vy4.*alpha));
S_vykappa = D_vykappa.*sin(r_vy5.*atan(r_vy6.*kappa)).*lambda_vykappa;
G_ykappa = cos(C_ykappa.*atan(B_ykappa.*kappa_s - E_ykappa.*(B_ykappa.*kappa_s - ...
    atan(B_ykappa.*kappa_s))))./cos(C_ykappa.*atan(B_ykappa.*S_hykappa - E_ykappa...
    .*(B_ykappa.*S_hykappa - atan(B_ykappa.*S_hykappa))));

%Lateral Force
F_y = D_y.*sin(C_y.*atan(B_y.*alpha_y-E_y.*(B_y.*alpha_y-atan(B_y.*alpha_y))))+S_vy;
F_y = transpose(G_ykappa.*F_y+S_vykappa);
Mz_Fycell = {S_hy,S_vy,K_yalpha,B_y,C_y,F_y,G_ykappa};

end
