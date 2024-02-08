%from NASA's interpolation on temp/pressure with altitude


X = 251460;                                                  %defining the altitude limit
rho = zeros(1,X);                                           %creating empty rho vector
T = zeros(1,X);                                          %temp at altitude - in C
P = zeros(1,X);

for z = 1:X 
    if z < 10999                                            %inside the troposphere
        T(z) = 15.04 - 0.00649 * z;                        
        P(z) = (101.29) * ((T(z) + 273.15)/288.08).^(5.256);
        rho(z) =  P(z)./(0.2869*(T(z) + 273.15));
    end
    if (11000 < z) && ( z < 24999)                          %inside the the lower stratosphere
        T(z) = -56.46;
        P(z) = 22.65 * exp(1.73 - 0.000157 * z);
        rho(z) = P(z)./(0.2869*(T(z) + 273.15));
    end
    if z > 25000                                            %upper stratosphere and beyond
        T(z) = -131.21 + 0.00299 * z;
        P(z) =  2.488 * ((T(z) + 273.15)/216.6).^(-11.388);
        rho(z) = P(z)./(0.2869*(T(z) + 273.15));
    end 
    rho = rho.';                                            %transpose of rho vector
    T = T.';
    P = P.';
    L = length(rho);
    z = 1:L;                                                %re-defining altitude vector to match rho vector dimensions                                                     
end
plot(z,rho);                                                %plotting rho as a function of z
Talt = T + 273*ones(1,L); %convert to Kelvin
Palt = P;%pressure (kPa)
save('rho_v_z.mat','rho','Talt', 'Palt')