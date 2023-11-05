function Cl = CL(M,A)

% the default approach - V-2's CL profile interpolation
mach = 0.2:0.2:5.4; aoa = [2 4 6 8 10];

CL_profile = [
    [0.3 0.63 1.05 1.49 2.05];
    [0.3 0.65 1.025 1.5 2.06];
    [0.3 0.62 1 1.46 2];
    [0.29 0.6 0.95 1.38 1.9];
    [0.31 0.65 1.01 1.42 1.92];
    [0.355 0.7 1.1 1.55 2.05];
    [0.35 0.69 1.06 1.49 1.97];
    [0.31 0.67 1 1.4 1.88];
    [0.3 0.62 0.96 1.34 1.77];
    [0.299 0.6 0.9 1.28 1.69];
    [0.298 0.58 0.88 1.2 1.6];
    [0.27 0.545 0.84 1.15 1.5];
    [0.25 0.51 0.8 1.09 1.41];
    [0.245 0.495 0.77 1.03 1.35];
    [0.23 0.48 0.73 0.99 1.29];
    [0.21 0.47 0.7 0.94 1.21];
    [0.205 0.45 0.69 0.9 1.18];
    [0.2 0.43 0.66 0.88 1.11];
    [0.199 0.405 0.62 0.85 1.09];
    [0.198 0.4 0.605 0.81 1.05];
    [0.197 0.399 0.599 0.8 1.01];
    [0.195 0.395 0.59 0.79 1];
    [0.193 0.39 0.57 0.78 0.98];
    [0.19 0.37 0.56 0.76 0.96];
    [0.185 0.37 0.55 0.73 0.92];
    [0.185 0.37 0.54 0.72 0.905];
    [0.185 0.36 0.52 0.705 0.9]
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
    grid_matrix = [CL_profile(squeeze_m_ind(1),squeeze_aoa_ind(1)) CL_profile(squeeze_m_ind(1),squeeze_aoa_ind(2));
        CL_profile(squeeze_m_ind(2),squeeze_aoa_ind(1)) CL_profile(squeeze_m_ind(2),squeeze_aoa_ind(2))];
    Cl = 1/diff(mach(squeeze_m_ind))/diff(aoa(squeeze_aoa_ind))*diff_mach*grid_matrix*diff_aoa;
elseif M > mach(end)
    Cl = 0;
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
    grid_matrix = [CL_profile(squeeze_m_ind(1),squeeze_aoa_ind(1)) CL_profile(squeeze_m_ind(1),squeeze_aoa_ind(2));
        CL_profile(squeeze_m_ind(2),squeeze_aoa_ind(1)) CL_profile(squeeze_m_ind(2),squeeze_aoa_ind(2))];
    Cl = 1/diff(mach(squeeze_m_ind))/diff(aoa(squeeze_aoa_ind))*diff_mach*grid_matrix*diff_aoa;
else
    Cl = 0;
end

% figure; plot(mach, CL_profile)
% [AOA,MACH] = meshgrid(aoa,mach);
% figure; surf(AOA,MACH,CL_profile); hold on;
% f1 = fit([AOA(:) MACH(:)],CL_profile(:),'poly45');
% machFit = 0.2:0.2:5.4; aoaFit = [0 4 6 8 10];
% [AOAFit, MachFit] = meshgrid(aoaFit,machFit);
% surf(AOAFit,MachFit,f1(AOAFit,MachFit));

end
