%% FUNCTION - single RK4 step
function [next_state,forces] = singleRK4(dynFunc,current_time,dt,current_state,current_control)
[k1,F1] = dynFunc(current_time,current_state,current_control);
[k2,F2] = dynFunc(current_time+dt/2,current_state+k1*dt/2,current_control);
[k3,F3] = dynFunc(current_time+dt/2,current_state+k2*dt/2,current_control);
[k4,F4] = dynFunc(current_time+dt,current_state+k3*dt,current_control);
next_state = current_state + 1/6*dt*(k1+2*k2+2*k3+k4);
forces = 1/6*(F1+2*F2+2*F3+F4);
