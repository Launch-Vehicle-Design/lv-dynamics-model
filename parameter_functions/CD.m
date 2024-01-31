function Cd = CD(M,A)

% the default approach - V-2's CD profile interpolation
mach = 0.2:0.2:5.4; aoa = [0 4 6 8 10];

CD_profile = [
    [0.15 0.18 0.24 0.3 0.4];
    [0.15 0.18 0.24 0.3 0.4];
    [0.15 0.18 0.24 0.3 0.4];
    [0.175 0.2 0.26 0.34 0.42];
    [0.3 0.35 0.4 0.5 0.65];
    [0.405 0.47 0.54 0.645 0.76];
    [0.35 0.4 0.475 0.575 0.7];
    [0.3 0.35 0.41 0.51 0.64];
    [0.26 0.3 0.375 0.46 0.58];
    [0.249 0.29 0.36 0.45 0.56];
    [0.24 0.28 0.35 0.44 0.55];
    [0.22 0.27 0.34 0.42 0.53];
    [0.21 0.26 0.33 0.4 0.5];
    [0.205 0.255 0.31 0.39 0.49];
    [0.2 0.25 0.3 0.38 0.475];
    [0.199 0.249 0.299 0.36 0.45];
    [0.19 0.24 0.28 0.35 0.44];
    [0.185 0.23 0.275 0.345 0.42];
    [0.18 0.215 0.26 0.325 0.4];
    [0.175 0.205 0.255 0.31 0.38];
    [0.17 0.2 0.25 0.3 0.37];
    [0.165 0.199 0.249 0.29 0.35];
    [0.16 0.19 0.235 0.275 0.34];
    [0.155 0.19 0.225 0.27 0.325];
    [0.15 0.175 0.21 0.255 0.305];
    [0.15 0.175 0.205 0.25 0.3];
    [0.15 0.175 0.2 0.249 0.299]
];

if M < mach(end) && M > mach(1) && A < aoa(end) && A > aoa(1)
    [~, adjac_m_ind] = min(abs(M-mach));
    if M-mach(adjac_m_ind) < 0
        squeeze_m_ind = [adjac_m_ind-1 adjac_m_ind];
    else 
        squeeze_m_ind = [adjac_m_ind adjac_m_ind+1];
    end
    [~, adjac_aoa_ind] = min(abs(A-aoa));
    if A-aoa(adjac_aoa_ind) < 0
        squeeze_aoa_ind = [adjac_aoa_ind-1 adjac_aoa_ind];
    else 
        squeeze_aoa_ind = [adjac_aoa_ind adjac_aoa_ind+1];
    end
    % bilinear interpolation
    m1 = mach(squeeze_m_ind(1)); m2 = mach(squeeze_m_ind(2));
    aoa1 = aoa(squeeze_aoa_ind(1)); aoa2 = aoa(squeeze_aoa_ind(2));
    diff_mach = [m2-M M-m1]; diff_aoa = [aoa2-A; A-aoa1];
    grid_matrix = [CD_profile(squeeze_m_ind(1),squeeze_aoa_ind(1)) CD_profile(squeeze_m_ind(1),squeeze_aoa_ind(2));
        CD_profile(squeeze_m_ind(2),squeeze_aoa_ind(1)) CD_profile(squeeze_m_ind(2),squeeze_aoa_ind(2))];
    Cd = 1/diff(mach(squeeze_m_ind))/diff(aoa(squeeze_aoa_ind))*diff_mach*grid_matrix*diff_aoa;
elseif M > mach(end)
    Cd = 0.2;
elseif A > aoa(end)
    [~, adjac_m_ind] = min(abs(M-mach));
    if M-mach(adjac_m_ind) < 0
        squeeze_m_ind = [adjac_m_ind-1 adjac_m_ind];
    else 
        squeeze_m_ind = [adjac_m_ind adjac_m_ind+1];
    end
    squeeze_aoa_ind = [length(aoa)-1 length(aoa)];
    % bilinear extrapolation
    m1 = mach(squeeze_m_ind(1)); m2 = mach(squeeze_m_ind(2));
    aoa1 = aoa(squeeze_aoa_ind(1)); aoa2 = aoa(squeeze_aoa_ind(2));
    diff_mach = [m2-M M-m1]; diff_aoa = [aoa2-A; A-aoa1];
    grid_matrix = [CD_profile(squeeze_m_ind(1),squeeze_aoa_ind(1)) CD_profile(squeeze_m_ind(1),squeeze_aoa_ind(2));
        CD_profile(squeeze_m_ind(2),squeeze_aoa_ind(1)) CD_profile(squeeze_m_ind(2),squeeze_aoa_ind(2))];
    Cd = 1/diff(mach(squeeze_m_ind))/diff(aoa(squeeze_aoa_ind))*diff_mach*grid_matrix*diff_aoa;
else 
    Cd = 0.2;
end

% figure; plot(mach, CD)
% figure; surf(AOA,MACH,CD); hold on;
% f1 = fit([AOA(:) MACH(:)],CD(:),'poly45');
% machFit = 0.2:0.2:5.4; aoaFit = [0 4 6 8 10];
% [AOAFit, MachFit] = meshgrid(aoaFit,machFit);
% surf(AOAFit,MachFit,f1(AOAFit,MachFit));

end
