function [history, perIndex] = dsouzaGuidLaw(initState,cond)
    % init history struct
    history.index = 1;
    history.time = nan(1e4,1);
    history.a = nan(1e4,3);
    history.v = nan(1e4,3);
    history.x = nan(1e4,3);
    
    time = 0; state = initState;
    while abs(state(1)) > cond.tolerance || abs(state(2)) > cond.tolerance || abs(state(3)) > cond.tolerance
        % control determination
        [t, control] = tgoSolve(state, cond);
        % history logging
        history.time(history.index) = time;
        history.a(history.index,:) = control';
        history.v(history.index,:) = state(4:6)';
        history.x(history.index,:) = state(1:3)';
        history.index = history.index + 1;
        
        % dynamical system update
        time = time + cond.dt;
        state = update(state,cond,control);
    end
    
    % figure; plot3(history.x(:,1),history.x(:,2),history.x(:,3));
    % performance index
    anorm = vecnorm(history.a')';
    anorm(all(isnan(anorm),2))=[];
    perIndex = 1/2*sum(anorm.^2*cond.dt);

    history.dv = sum(anorm*cond.dt);
end
    
%% Inverse Square Dynamics Updates
function stateNext = update(state, cond, av)
    xv = state(1:3); uv = state(4:6);
    a = av + [0;0;(-cond.mu/(xv(3)+cond.R)^2)];
    xvNext = xv + uv.*cond.dt + 0.5*a.*cond.dt^2;
    uvNext = uv + a.*cond.dt;
    stateNext = [xvNext; uvNext];
end

%% Time-to-Go Solver and Control
function [tgo, control] = tgoSolve(state, cond)
    r = state(1:3); v = state(4:6);
    a = -2*dot(v,v)/(cond.gamma+(-cond.mu/(cond.R)^2)^2/2);
    b = -12*dot(v,r)/(cond.gamma+(-cond.mu/(cond.R)^2)^2/2);
    c = -18*dot(r,r)/(cond.gamma+(-cond.mu/(cond.R)^2)^2/2);
    
    alpha = 1/3*(3*(a^2-4*c)-4*a^2);
    beta = 1/27*(16*a^3-18*a*(a^2-4*c)-27*b^2);
    delta = alpha^3/27+beta^2/4;
    if delta > 0
        Z = nthroot(-beta/2+sqrt(delta),3) + nthroot(-beta/2-sqrt(delta),3);
    else
        Z = (-beta/2+sqrt(delta))^(1/3)+(-beta/2-sqrt(delta))^(1/3);
    end
    eta = Z - 2*a/3;
    zeta = -b/2/eta; xi = (a+eta)/2;
    
    tgo1 = (sqrt(eta) + sqrt(eta-4*(xi-sqrt(eta)*zeta)))/2;
    tgo2 = (sqrt(eta) - sqrt(eta-4*(xi-sqrt(eta)*zeta)))/2;
    tgo3 = (-1*sqrt(eta) + sqrt(eta-4*(xi+sqrt(eta)*zeta)))/2;
    tgo4 = (-1*sqrt(eta) - sqrt(eta-4*(xi+sqrt(eta)*zeta)))/2;
    if (imag(tgo1)~=0 || tgo1<0) 
        tgo1 = 0; end
    if (imag(tgo2)~=0 || tgo2<0) 
        tgo2 = 0; end
    if (imag(tgo3)~=0 || tgo3<0) 
        tgo3 = 0; end
    if (imag(tgo4)~=0 || tgo4<0) 
        tgo4 = 0; end
    
    % tgo_fzero = fzero(@(tgo) tgo^4+a*tgo^2+b*tgo+c,10);
    tgo = tgo1 + tgo2 + tgo3 + tgo4;
    control = -4*v/tgo-6*r/tgo^2-[0 0 (-cond.mu/(cond.R)^2)]';
end