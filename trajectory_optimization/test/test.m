clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

file_name = "dircol_sqp_final1_scaled.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end
states = nan([4,500]);
t = 1:500;
for i = 1:500 states(:,i) = continuous_state_interpreter(i,log_x,log_param); end

subplot(2,2,1); plot(t, states(1,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,2,2); plot(t, states(2,:)*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
subplot(2,2,3); plot(t, states(3,:)/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Altitude (km)");
subplot(2,2,4); plot(t, states(4,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");

velo = states(1,:); alti = states(3,:);
q = (1/2*log_param.rho(alti./log_param.scales.length).*(velo./log_param.scales.speed).^2).*log_param.scales.pressure;
figure; plot(t,q,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Dynamics Pressure (Pa)");

save("alti.mat","t","alti");