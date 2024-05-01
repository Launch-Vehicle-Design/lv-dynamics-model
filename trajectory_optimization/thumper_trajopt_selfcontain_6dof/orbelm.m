function elements = orbelm(r,v,param)
%% 6 orbit elements
i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

mu = param.mu*param.scales.gravparam;

% semi major
specificEnergy = norm(v)^2 / 2 - mu / norm(r);
a = -1*mu / (2*specificEnergy);

% eccentricity
h = cross(r,v);
e = cross(v,h) / mu - r/norm(r);
eccentricity = norm(e);

% inclination
inclination = acos(dot(k,h) / norm(h)) * 180/pi;

% right ascension of the ascending node
n = cross(k,h) / norm(cross(k,h)); % unit vector in line of node
raan = acos(dot(n,i)) * 180/pi;
if dot(n,j) < 0
    raan = 360 - raan;
end

% argument of periapsis
omega = acos(dot(e,n) / norm(e)) * 180/pi;
if dot(e,k) < 0
    omega = 360 - omega;
end

% true anomaly
nu = acos(dot(e,r) / (norm(e)*norm(r))) * 180/pi;
if dot(r,v) < 0
    nu = 360 - nu;
end

elements.a = a;
elements.e_vec = e;
elements.e = eccentricity;
elements.inc = inclination;
elements.raan = raan;
elements.omega = omega;
elements.nu = nu;
elements.mu = mu;
end

