function [history, perIndex] = linearGuidLaw(initState,cond)
    % init history struct
    history.index = 1;
    history.time = nan(1e4,1);
    history.a = nan(1e4,3);
    history.v = nan(1e4,3);
    history.x = nan(1e4,3);

    time = 0;
    % convert to Local Surface XYZ
    state = DCA2LSXYZ(initState,cond);
    [lambdai, lambdavi] = lambdas();
    while abs(state(1)) > cond.tolerance || abs(state(2)) > cond.tolerance || abs(state(3)) > cond.tolerance
        % control determination
        [t, control] = tgoSolve(state, cond, lambdai, lambdavi);
        cond.tgoGuess = t;

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
        
    history.stateDCA = LSXYZ2DCA([history.x history.v], cond);
    % figure; plot3(history.stateDCA(:,1),history.stateDCA(:,2),history.stateDCA(:,3));
    % performance index
    anorm = vecnorm(history.a')';
    anorm(all(isnan(anorm),2))=[];
    perIndex = 1/2*sum(anorm.^2*cond.dt);
    
    history.dv = sum(anorm*cond.dt);
end

%% Inverse Square Dynamics Updates
function stateNext = update(state, cond, av)
    xv = state(1:3) + cond.targetLocation; uv = state(4:6);
    a = -cond.mu/(xv'*xv)^1.5*xv + av;
    xvNext = xv + uv.*cond.dt + 0.5*a.*cond.dt^2;
    uvNext = uv + a.*cond.dt;
    stateNext = [xvNext - cond.targetLocation; uvNext];
end

%% Downrange-Crossrange-Altitude to Local Surface XYZ
function stateLSXYZ = DCA2LSXYZ(state, cond)
    x = state(1:3); v = state(4:6);
    % position conversion - DCA -> LocalSpherical -> MCI -> LSXYZ
    theta = x(1)/cond.R; phi = x(2)/cond.R; r = cond.R + x(3);
    xMCI = [r*cos(phi)*cos(theta); r*cos(phi)*sin(theta); r*sin(phi)];
    xLSXYZ = xMCI - cond.targetLocation;
    
    % velocity conversion - rotation matrix using LocalSpherical
    C = [-sin(theta) cos(theta) 0;...
        -sin(phi)*cos(theta) -sin(phi)*sin(theta) cos(phi);...
        cos(phi)*cos(theta) cos(phi)*sin(theta) sin(phi)]';
    vLSXYZ = C*v;
    stateLSXYZ = [xLSXYZ; vLSXYZ];
end

%% Time-to-Go Solve Function Definition
function [lambdai, lambdavi] = lambdas()
    % construct functions
    cti = @(tgo,Gi) cosh(sqrt(Gi)*tgo); sti = @(tgo,Gi) sinh(sqrt(Gi)*tgo);
    Ci = @(tgo,Gi,Gci) Gci/Gi*(cti(tgo, Gi)-1); Cvi = @(tgo,Gi,Gci) -Gci/sqrt(Gi)*sti(tgo, Gi);
    d1 = @(tgo,Gi,Gci) sti(tgo,Gi)/2/Gi^1.5-tgo*cti(tgo,Gi)/2/Gi; d2 = @(tgo,Gi,Gci) -tgo*sti(tgo,Gi)/2/sqrt(Gi);
    d3 = @(tgo,Gi,Gci) -d2(tgo,Gi,Gci); d4 = @(tgo,Gi,Gci) tgo*cti(tgo, Gi)/2+sti(tgo, Gi)/2/sqrt(Gi);
    nui = @(tgo,Gi,Gci,x,u) (d4(tgo,Gi,Gci)*(x-Ci(tgo,Gi,Gci))-d2(tgo,Gi,Gci)*(u-Cvi(tgo,Gi,Gci)))...
        /(d1(tgo,Gi,Gci)*d4(tgo,Gi,Gci)-d2(tgo,Gi,Gci)*d3(tgo,Gi,Gci));
    nuvi = @(tgo,Gi,Gci,x,u) (-d3(tgo,Gi,Gci)*(x-Ci(tgo,Gi,Gci))+d1(tgo,Gi,Gci)*(u-Cvi(tgo,Gi,Gci)))...
        /(d1(tgo,Gi,Gci)*d4(tgo,Gi,Gci)-d2(tgo,Gi,Gci)*d3(tgo,Gi,Gci));
    lambdai = @(tgo,Gi,Gci,x,u) 1/2*(nui(tgo,Gi,Gci,x,u)-sqrt(Gi)*nuvi(tgo,Gi,Gci,x,u))*exp(-sqrt(Gi)*tgo)...
        + 1/2*(nui(tgo,Gi,Gci,x,u)+sqrt(Gi)*nuvi(tgo,Gi,Gci,x,u))*exp(sqrt(Gi)*tgo);
    lambdavi = @(tgo,Gi,Gci,x,u) -1/2/sqrt(Gi)*(nui(tgo,Gi,Gci,x,u)-sqrt(Gi)*nuvi(tgo,Gi,Gci,x,u))*exp(-sqrt(Gi)*tgo)...
        + 1/2/sqrt(Gi)*(nui(tgo,Gi,Gci,x,u)+sqrt(Gi)*nuvi(tgo,Gi,Gci,x,u))*exp(sqrt(Gi)*tgo);
end

%% Time-to-Go Solver and Control
function [t, control] = tgoSolve(state, cond, lambdai, lambdavi)
    sum = @(tgo) 0; a = @(tgo) [0; 0; 0];
    for i = 1:3
        Gi = cond.Gi(i,i); Gci = cond.Gc(i); 
        x = state(i); u = state(i+3);
        sum = @(tgo) sum(tgo) - lambdavi(tgo,Gi,Gci,x,u)^2/2+lambdai(tgo,Gi,Gci,x,u)*u...
            +lambdavi(tgo,Gi,Gci,x,u)*(Gci+Gi*x);
        switch i
            case 1
                a = @(tgo) a(tgo) + [-1*lambdavi(tgo,Gi,Gci,x,u); 0; 0];
            case 2
                a = @(tgo) a(tgo) + [0; -1*lambdavi(tgo,Gi,Gci,x,u); 0];
            case 3
                a = @(tgo) a(tgo) + [0; 0; -1*lambdavi(tgo,Gi,Gci,x,u)];
        end
    end
    t = fzero(sum,cond.tgoGuess,cond.options);
    control = a(t);
end

%% Unwrap the Result from Target Center Inertial to Downrange-Crossrange-Altitude
function stateDCA = LSXYZ2DCA(state, cond)
    x = state(:,1:3); u = state(:,4:6);
    % position conversion - LSXYZ/TCI -> MCI/CCI -> LocalSpherical -> DCA
    % celestial body center inertial
    xCCI = x + repmat(cond.targetLocation',size(x,1),1);
    rXY = vecnorm(xCCI(:,1:2)')'; r = vecnorm(xCCI')';
    phi = asin(xCCI(:,3)./rXY); theta = atan2(xCCI(:,2), xCCI(:,1));
    xTCI = [cond.R.*theta cond.R.*phi r-cond.R];
    
    % velocity conversion - rotation matrix using LocalSpherical
    
    stateDCA = [xTCI];
end
