clear;
close all;
fclose all;
clc;

% PART 1
mdot = 75; %kg/s
Mz2 = 0.5;
p02 = 35000; %Pa
T02 = 250; %K
pi_c = 25;
Mt = 1.4;
ec = 0.95;
gamma = 1.4;
R = 287; %J/kgK
zeta_f = 0.4:0.01:0.6;
cp = 1000; %J/kgK

% a) Calculate Af
Af = mdot / (Mz2 * p02) * sqrt((R * T02) / gamma) * ...
    (1 + (gamma - 1) / 2 * Mz2^2)^((gamma + 1)/(2 * (gamma - 1)));
fprintf('The compressor face flow-through area, A_f, is %f m^2\n', Af);

% b) Find Blade tip radius & blade length vs. zeta_f
for i = 1:size(zeta_f, 2)
    r_tf(i) = sqrt(Af / (pi * (1 - zeta_f(i)^2)));
    l_f(i) = r_tf(i) * (1 - zeta_f(i));
end
figure(1);
plot(zeta_f, r_tf, zeta_f, l_f,'--');
xlabel('\zeta_f');
ylabel('Length [m]');
legend('blade tip radius: r_{t,f}', 'blade length: l_f','Location','best')
title('Blade tip radius & Blade length vs. \zeta_f');
xlim([0.4 0.6]);

% c) Find RPM & U_pitch vs. zeta_f
for i=1:size(zeta_f, 2)
    omega(i) = Mt * sqrt(gamma * R * T02 / ...
        ((1 + (gamma - 1) / 2 * Mz2^2))) / r_tf(i);
    RPM(i) = omega(i) / (2 * pi) * 60;
    U_pf(i) = omega(i) * r_tf(i) * (1 + zeta_f(i)) / 2;
end
figure(2);
plot(zeta_f, U_pf);
xlabel('\zeta_f');
ylabel('U_{p,f} [m/s]');
title('U_{pf} vs. \zeta_f');
xlim([0.4 0.6]);

figure(3);
plot(zeta_f, RPM);
xlabel('\zeta_f');
ylabel('\Omega [RPM]');
title('\Omega vs. \zeta_f');
xlim([0.4 0.6]);

% At zeta_f = 0.5, the values are:
U_pf_ = U_pf(11);
omega_ = omega(11);
rpm_ = RPM(11);
r_tf_ = r_tf(11);
l_f_ = l_f (11);

fprintf('\nAt zeta = 0.5\n--------------------------------\n');
fprintf('U_{pf} = %f m/s\nOmega = %f rad/s OR %f RPM\nr_{t,f} = %f m\nl_f = %f m\n',...
    U_pf_, omega_, rpm_, r_tf_, l_f_);


% PART 2
% Some design parameters:
sigma_r = 1.0; % rotor solidity
sigma_s = 1.25; % stator solidity
cz = Mz2 * sqrt(gamma * R * T02 / (1 + (gamma - 1) / 2 * Mz2^2)); % const throughout
c3 = cz; % const throughout
zeta_f = 0.5; % for first stage only

fprintf('C_z = %f m/s\n\n', cz);

% a)
% Degree of reaction range - pick one
DOR = 0.6:0.00001:1.0;

% First stage
Mz(1) = Mz2;
T0(1) = T02;
P0(1) = p02;
U(1) = Mt * (1 + zeta_f)/2 * sqrt(gamma * R * T02 / ...
        (1 + (gamma - 1) / 2 * Mt^2)); %m/s
U(1) = -U(1); % direction is negative in velocity triangle
Af(1) = Af;
r_t(1) = r_tf_;
r_h = zeta_f * r_t(1);  % CONSTANT FOR THIS PROBLEM

% BONUS:
%r_h(1) = zeta_f * r_t(1);
%r_p  = (r_t(1) + r_h(1))/2          % CONSTANT FOR THIS PROBLEM

om = omega_; % constant


% convention:
% beta, wtheta +
% alpha, U -
beta1(1) = atan(-U(1) / cz);
wtheta1(1) = cz * tan(beta1(1));
w1(1) = sqrt(cz^2 + wtheta1(1)^2);

% Method 2: Another way of defining velocity triangles (checking 2 ways for same ans):
% wtheta1(1) = -U(1);
% w1(1) = sqrt(cz^2 + wtheta1(1)^2);
 alpha1 = 0; % const
% beta1(1) = atan(wtheta1(1) / cz);

for i = 1:size(DOR, 2)
    % Defining first stage with acceptable Dfactor

    beta2 = atan(2 * U(1) * (0.5 - DOR(i)) / cz); % const
    wtheta2(1) = cz * tan(beta2);
    w2(1) = sqrt(cz^2 + wtheta2(1)^2);
    alpha2(1) = atan((U(1) + wtheta2(1)) / cz);
    dctheta(1) = cz * tan(alpha2(1));
    c2(1) = sqrt(cz^2 + dctheta(1)^2);

% Method 2: Testing different ways of calculating velocity triangles:
%     wtheta2(1) = -2 * DOR(i) * U(1) - wtheta1(1);
%     beta2 = atan(wtheta2(1) / cz); % const
%     dctheta(1) = U(1) + cz * tan(beta2); % = wtheta2 - wtheta1
%     w2(1) = sqrt(cz^2 + wtheta2(1)^2);
%     c2(1) = sqrt(cz^2 + (dctheta(1))^2);
%     alpha2(1) = atan((U(1) + cz * tan(beta2)) / cz); % first stage
    
    % Calculate Dfactor
    Dr(1) = 1 - w2(1) / w1(1) + abs(wtheta1(1) - wtheta2(1)) / (2 * sigma_r * w1(1));
    Ds(1) = 1 - c3 / c2(1) + abs(dctheta(1)) / (2 * sigma_s * c2(1));
    done = 0;
    if (Dr(1) < 0.6 && Ds(1) < 0.6)
        DORmin = DOR(i); % choose a DOR that satisfy Dfactor
        
        % for stage 2 and beyond
        j = 2;
        % clear arrays after 2nd element from before
        T0(2:size(T0,2)) = []; 
        P0(2:size(P0,2)) = [];
        pic = [];
        Mz(2:size(Mz, 2)) = [];
        Af(2:size(Af,2)) = [];
        r_t(2:size(r_t,2)) = [];
        om(2:size(om,2)) = [];
        U(2:size(U,2)) = [];
        wtheta1(2:size(wtheta1,2)) = [];
        w1(2:size(w1,2)) = [];
        beta1(2:size(beta1,2)) = [];
        wtheta2(2:size(wtheta2,2)) = [];
        alpha2(2:size(alpha2,2)) = [];
        dctheta(2:size(dctheta,2)) = [];
        w2(2:size(w2,2)) = [];
        c2(2:size(c2,2)) = [];
        Dr(2:size(Dr,2)) = [];
        Ds(2:size(Ds,2)) = [];
        picsofar(1) = 1;
        
        while (1) % while pic is less than 25
            T0(j) = T0(j - 1) * (1 + U(j - 1) * dctheta(j - 1)/(cp * T0(j - 1)));
            T0ratio(j-1) = (1 + U(j - 1) * dctheta(j - 1)/(cp * T0(j - 1)));
            P0(j) = P0(j - 1) * (T0(j) / T0(j - 1))^(gamma * ec / (gamma - 1));
            pic(j - 1) = (T0(j) / T0(j - 1))^(gamma * ec / (gamma - 1));
            picsofar = picsofar * pic(j-1);

            
            Mz(j) = cz / sqrt(gamma * R * T0(j) - (gamma - 1)/2 * cz^2);
            Af(j) = mdot / (Mz(j) * P0(j)) * sqrt((R * T0(j)) / gamma) * ...
                 (1 + (gamma - 1) / 2 * Mz(j)^2)^((gamma + 1)/(2 * (gamma - 1)));
            r_t(j) = sqrt(Af(j)/pi + r_h^2);
            U(j) = om * (r_t(j) + r_h) / 2;
            U(j) = -U(j);
            
            % calculate the velocity triangles

            beta1(j) = atan(-U(j)/cz);
            wtheta1(j) = cz*tan(beta1(j));
            wtheta2(j) = cz*tan(beta2);
            w1(j) = sqrt(cz^2 + wtheta1(j)^2);
            w2(j) = sqrt(cz^2 + wtheta2(j)^2);
            alpha2(j) = atan((U(j) + wtheta2(j))/cz);
            dctheta(j) = cz*tan(alpha2(j));
            c2(j) = sqrt(cz^2 + dctheta(j)^2);
    
% Method 2:            
%             wtheta1(j) = -U(j);
%             w1(j) = sqrt(cz^2 + wtheta1(j)^2);
%             beta1(j) = atan(wtheta1(j) / cz);
%             wtheta2(j) = cz * tan(beta2);
%             alpha2(j) = atan((U(j) + cz * tan(beta2)) / cz);
%             dctheta(j) = U(j) + cz * tan(beta2);
%             w2(j) = sqrt(cz^2 + wtheta2(j)^2);
%             c2(j) = sqrt(cz^2 + (dctheta(j))^2);
            
            % At every stage check Dfactor
            Dr(j) = 1 - w2(j) / w1(j) + abs(wtheta1(j) - wtheta2(j)) /...
                (2 * sigma_r * w1(j));
            Ds(j) = 1 - c3 / c2(j) + abs(dctheta(j)) / (2 * sigma_s * c2(j));
            if (Dr(j) >= 0.6 || Ds(j) >= 0.6)
                break; % choose another DOR because this one stalls!
            end
            
            if (picsofar > pi_c) % higher than 25, finish!
                done = 1;
                break;
            end
            
            j = j + 1;
        end
            
      if (done) % process is done, no need to look for another DOR
         break;
      end
    end
end
if(~done)
    disp ('Sorry! Could not find a DOR!')
    return;
end

% a) DOR
fprintf('Minimum degree of reaction needed to avoid stall = %.2f\n', DORmin);
% b) Number of stages and compression ratio achieved
fprintf('Number of stages = %d\nCompression Ratio achieved = %.2f\n',j - 1, picsofar);

% c) Velocity triangles
% First stage
fprintf('\nFirst Stage\n------------------------\n');
fprintf('wtheta1 = %f m/s\nwtheta2 = %f m/s\ndctheta = %f m/s\nw1 = %f m/s\nw2 = %f m/s\n', ...
    wtheta1(1),wtheta2(1),abs(dctheta(1)),w1(1),w2(1));
fprintf('c1 = cz = c3 = %f m/s\nc2 = %f m/s\nU = %f m/s\nalpha1 = alpha3 = 0 deg\nalpha2 = %f deg\n', ...
    cz, c2(1),U(1),alpha2(1)*180/pi);
fprintf('beta1 = %f deg\nbeta2 = %f deg\n', ...
    beta1(1)*180/pi,beta2*180/pi);

fprintf('\nMiddle Stage (Stage 8)\n------------------------\n');
fprintf('wtheta1 = %f m/s\nwtheta2 = %f m/s\ndctheta = %f m/s\nw1 = %f m/s\nw2 = %f m/s\n', ...
    wtheta1(8),wtheta2(8),abs(dctheta(8)),w1(8),w2(8));
fprintf('c1 = cz = c3 = %f m/s\nc2 = %f m/s\nU = %f m/s\nalpha1 = alpha3 = 0 deg\nalpha2 = %f deg\n', ...
    cz, c2(8),U(8),alpha2(8)*180/pi);
fprintf('beta1 = %f deg\nbeta2 = %f deg\n', ...
    beta1(8)*180/pi,beta2*180/pi);

fprintf('\nLast Stage (Stage 15)\n------------------------\n');
fprintf('wtheta1 = %f m/s\nwtheta2 = %f m/s\ndctheta = %f m/s\nw1 = %f m/s\nw2 = %f m/s\n', ...
    wtheta1(16),wtheta2(16),abs(dctheta(16)),w1(16),w2(16));
fprintf('c1 = cz = c3 = %f m/s\nc2 = %f m/s\nU = %f m/s\nalpha1 = alpha3 = 0 deg\nalpha2 = %f deg\n', ...
    cz, c2(16),U(16),alpha2(16)*180/pi);
fprintf('beta1 = %f deg\nbeta2 = %f deg\n', ...
    beta1(16)*180/pi,beta2*180/pi);

% d) Compressor conditions
% Stagnation pressure
figure(4);
plot(0:j-1,P0);
title('Stagnation Pressure');
ylabel('p_0 [Pa]');
xlabel('Stage Number');
% Stagnation temperature
figure(5);
plot(0:j-1,T0);
title('Stagnation Temperature');
ylabel('T_0 [K]');
xlabel('Stage Number');
% Stagnation pressure ratio
figure(6);
plot(1:j-1,pic);
title('Stagnation Pressure Ratio');
ylabel('P_0_i/P_0_{i-1} [Pa]');
xlabel('Stage Number');
xlim([1 15]);
% Stagnation temperature ratio
figure(7);
plot(1:j-1,T0ratio);
title('Stagnation Temperature Ratio');
ylabel('T_0_i/T_0_{i-1} [K]');
xlabel('Stage Number');
xlim([1 15]);
% Dfactor
figure(8);
plot(0:j-1,Dr, 0:j-1, Ds, '--');
title('D-factor');
ylabel('D-factor');
xlabel('Stage Number');
legend('D_r','D_s');

% e)
%Blade Tip
%Blade Hub
%Pitch-line radius
r_h(1:j) = r_h;
figure(9);
plot(0:j-1,r_t, 0:j-1,r_h, 0:j-1,(r_t+r_h)/2, '--');
title('Radius of compressor parameters');
ylabel('Radius [m]');
xlabel('Stage Number');
legend('Blade tip', 'Blade Hub', 'Pitch-line Radius');




