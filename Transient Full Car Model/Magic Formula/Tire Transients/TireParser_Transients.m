function [SA_out, FY_out, FZ_out, MZ_out, MX_out, IA_out, P_out, index2, Time_out, V_out]...
    = TireParser_Transients(P_input, IA_input, FZ_input, V_input, SA_input)

% P_input = tire pressure
% IA_input = inclination angle (camber) 
% FZ_input = normal force in z direction 
% V_input = velocity in mph
% SA_input = slip angle in degrees

P_tol = .8;
IA_tol = 0.5;
FZ_tol = 60;
%V_tol = 25;
SA_tol = 0.9;

FZ_input = -FZ_input; %convention

load('A1654run32.mat');
%load('A1965run14.mat');
%load('A1965run5.mat');

index = [];
index2 = [];
index1 = [];
index3 = [];
%isolate data based on given parameters

for i = 1:numel(P_input)
    for j = 1:numel(IA_input)
        for k = 1:numel(FZ_input)
                for m = 1:numel(SA_input)
                    index2 = [index2; find(P > P_input(i) - P_tol & P < P_input(i) + P_tol & IA > IA_input(j) - IA_tol...
                        & IA < IA_input(j) + IA_tol & FZ > FZ_input(k) - FZ_tol & FZ < FZ_input(k) + FZ_tol & SA > SA_input(m) - SA_tol & SA < SA_input(m) + SA_tol &...
                        V > 0.11)];
                end
        end
    end
end

index2 = index2(1:end);
indexnew = index2(1:floor(numel(index2)*0.5));
%outputs
% SA_out = slip angle 
% FY_out = lateral force 
% FZ_out = normal force
% MZ_out = aligning torque
% MX_out = overturning moment 
% IA_out = inclination angle (camber) 
% P_out = tire pressure

SA_out = SA(indexnew);
FY_out = FY(indexnew);
FZ_out = FZ(indexnew);
MZ_out = MZ(indexnew);
MX_out = MX(indexnew);
IA_out = IA(indexnew);
P_out = P(indexnew);
Time_out = ET(indexnew);
Time_out = Time_out - Time_out(1);
FY_out = FY_out - FY_out(1);
V_out = V(indexnew);
end