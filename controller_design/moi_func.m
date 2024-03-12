%% FUNCTION - generate the Moment of Inertia calculation for each component
function funcs = moi_func
% variable naming explanation
% moi - moment of inertia
% c - (wrt) component center of mass
% t - (wrt) component top center line
% lc - (wrt) sub-component center of mass
% a - additional moment of inertia from parallel axis theorem

% parallel axis theorem
moia = @(m,r) m*(r'*r*eye(3) - r*r');
funcs.moia = moia;

% % % fixed structure % % %
% payload - cube
funcs.rc_pl = @(a) [-1/2*a; 0; 0];
funcs.moic_pl = @(mpl,a) mpl*a^2/6*eye(3);
funcs.moit_pl = @(mpl,a) funcs.moic_pl(mpl,a) + moia(mpl,-funcs.rc_pl(a));
% payload fairing - cone
funcs.rc_plf = @(d,l) [-2/3*l; 0; 0];
funcs.moit_plf = @(mplf,d,l) mplf*diag([1/8*d^2 1/16*d^2+1/6*l^2 1/16*d^2+1/6*l^2]);
funcs.moic_plf = @(mplf,d,l) funcs.moit_plf(mplf,d,l) - moia(mplf,-funcs.rc_plf(d,l));
% second stage propellant tanks structure - two stacked hollow cylinder with three domes
funcs.rc_2ndstg_td = @(d,lo,lf) [0; 0; 0];
funcs.rc_2ndstg_md = @(d,lo,lf) [-lf; 0; 0];
funcs.rc_2ndstg_ld = @(d,lo,lf) [-(lf+lo); 0; 0];
funcs.rc_2ndstg_ft = @(d,lo,lf) [-lf/2; 0; 0];
funcs.rc_2ndstg_ot = @(d,lo,lf) [-lf-lo/2; 0; 0];
funcs.rc_2ndstg_pt = @(d,lo,lf) [-(d*(lo+2*lf)+2*(lo+lf)^2)/(3*d+4*(lo+lf)); 0; 0];
funcs.moilc_dome = @(mps_dome,d) mps_dome*diag([1/8*d^2 1/16*d^2 1/16*d^2]);
funcs.moilc_tank = @(mps_tank,d,l) mps_tank*diag([1/4*d^2 1/8*d^2+1/12*l^2 1/8*d^2+1/12*l^2]);
funcs.rho_2ndstg_pt = @(mps,d,lo,lf) mps/(3/4*pi*d^2+pi*d*(lf+lo));
funcs.moit_2ndstg_pt = @(mps,d,lo,lf) funcs.rho_2ndstg_pt(mps,d,lo,lf) * ...
    (3*funcs.moilc_dome(1/4*pi*d^2,d) + ...
    moia(1/4*pi*d^2,-funcs.rc_2ndstg_td(d,lo,lf)) + ...
    moia(1/4*pi*d^2,-funcs.rc_2ndstg_md(d,lo,lf)) + ...
    moia(1/4*pi*d^2,-funcs.rc_2ndstg_ld(d,lo,lf)) + ...
    funcs.moilc_tank(pi*d*lf,d,lf) + funcs.moilc_tank(pi*d*lo,d,lo) + ...
    moia(pi*d*lf,-funcs.rc_2ndstg_ft(d,lo,lf)) + ...
    moia(pi*d*lo,-funcs.rc_2ndstg_ot(d,lo,lf)));
funcs.moic_2ndstg_pt = @(mps,d,lo,lf) funcs.moit_2ndstg_pt(mps,d,lo,lf) - moia(mps,-funcs.rc_2ndstg_pt(d,lo,lf));
% second stage engine structure - cone
funcs.rc_2ndstg_eng = @(d,l) [-2/3*l; 0; 0];
funcs.moit_2ndstg_eng = @(m_rdre,d,l) m_rdre*diag([1/8*d^2 1/16*d^2+1/6*l^2 1/16*d^2+1/6*l^2]);
funcs.moic_2ndstg_eng = @(m_rdre,d,l) funcs.moit_2ndstg_eng(m_rdre,d,l) - moia(m_rdre,-funcs.rc_2ndstg_eng(d,l));
% second stage GNC structure - rod for external wiring harness
funcs.rc_2ndgnc = @(d,l) [-l/2; 0; 0];
funcs.rlc_2ndgnc = @(d,l,angle) [-l/2; d/2*cos(angle); d/2*sin(angle)];
funcs.moit_2ndgnc = @(m_gnc,d,l) m_gnc*diag([0; 1/12*l^2; 1/12*l^2]) + ...
    moia(m_gnc,-funcs.rlc_2ndgnc(d,l,15/180*pi)) + moia(m_gnc,-funcs.rlc_2ndgnc(d,l,195/180*pi));
funcs.moic_2ndgnc = @(m_gnc,d,l) funcs.moit_2ndgnc(m_gnc,d,l) - moia(m_gnc,-funcs.rc_2ndgnc(d,l));
% interstage skirt - hollow cylinder
funcs.rc_intstg = @(d,l) [-l/2; 0; 0];
funcs.moic_intstg = @(m_intstg,d,l) m_intstg*diag([1/4*d^2 1/8*d^2+1/12*l^2 1/8*d^2+1/12*l^2]);
funcs.moit_intstg = @(m_intstg,d,l) funcs.moic_intstg(m_intstg,d,l) + moia(m_intstg,-funcs.rc_intstg(d,l));
% first stage solid casing - hollow cylinder with two domes
funcs.rc_1ststg_td = @(d,l) [0; 0; 0];
funcs.rc_1ststg_ld = @(d,l) [-l; 0; 0];
funcs.rc_1ststg_sc = @(d,l) [-l/2; 0; 0];
funcs.rc_1ststg_pt = @(d,l) [-(d*l+2*l^2)/(2*d+4*l); 0; 0];
funcs.rho_1ststg_pt = @(mps,d,l) mps/(pi*d^2/2+pi*d*l);
funcs.moit_1ststg_pt = @(mps,d,l) funcs.rho_1ststg_pt(mps,d,l) * ...
    (2*funcs.moilc_dome(1/4*pi*d^2,d) + ...
    moia(1/4*pi*d^2,-funcs.rc_1ststg_td(d,l)) + ...
    moia(1/4*pi*d^2,-funcs.rc_1ststg_ld(d,l)) + ...
    funcs.moilc_tank(pi*d*l,d,l) + ...
    moia(pi*d*l,-funcs.rc_1ststg_sc(d,l)));
funcs.moic_1ststg_pt = @(mps,d,l) funcs.moit_1ststg_pt(mps,d,l) - moia(mps,-funcs.rc_1ststg_pt(d,l));
% first stage engine structure - cone
funcs.rc_1ststg_eng = @(d,l) [-2/3*l; 0; 0];
funcs.moit_1ststg_eng = @(m_sm,d,l) m_sm*diag([1/8*d^2 1/16*d^2+1/6*l^2 1/16*d^2+1/6*l^2]);
funcs.moic_1ststg_eng = @(m_sm,d,l) funcs.moit_1ststg_eng(m_sm,d,l) - moia(m_sm,-funcs.rc_1ststg_eng(d,l));
% first stage GNC structure - rod for external wiring harness
funcs.rc_1stgnc = @(d,l) [-l/2; 0; 0];
funcs.rlc_1stgnc = @(d,l,angle) [-l/2; d/2*cos(angle); d/2*sin(angle)];
funcs.moit_1stgnc = @(m_gnc,d,l) m_gnc*diag([0; 1/12*l^2; 1/12*l^2]) + ...
    moia(m_gnc,-funcs.rlc_1stgnc(d,l,15/180*pi)) + moia(m_gnc,-funcs.rlc_1stgnc(d,l,195/180*pi));
funcs.moic_1stgnc = @(m_gnc,d,l) funcs.moit_1stgnc(m_gnc,d,l) - moia(m_gnc,-funcs.rc_1stgnc(d,l));
% aft skirt - hollow cylinder with one dome (low side)
funcs.rc_aft_d = @(d,l) [-l; 0; 0];
funcs.rc_aft_t = @(d,l) [-l/2; 0; 0];
funcs.rc_aft = @(d,l) [-(d*l+2*l^2)/(d+4*l); 0; 0];
funcs.rho_aft = @(m_aft,d,l) m_aft/(pi*d^2/4+pi*d*l);
funcs.moit_aft = @(m_aft,d,l) funcs.rho_aft(m_aft,d,l) * ...
    (funcs.moilc_dome(1/4*pi*d^2,d) + moia(1/4*pi*d^2,-funcs.rc_aft_d(d,l)) + ...
    funcs.moilc_tank(pi*d*l,d,l) + moia(pi*d*l,-funcs.rc_aft_t(d,l)));
funcs.moic_aft = @(m_aft,d,l) funcs.moit_aft(m_aft,d,l) - moia(m_aft,-funcs.rc_aft(d,l));

% % % propellant % % %
% second stage liquid propellant - spherical pendulum model
funcs.rc_2ndstg_p = @(d,l) [-l/2; 0; 0];
funcs.moic_2ndstg_p_rgb = @(mp,d,l) mp*diag([1/8*d^2 1/16*d^2+1/12*l^2 1/16*d^2+1/12*l^2]);
funcs.moit_2ndstg_p_rgb = @(mp,d,l) funcs.moic_2ndstg_p_rgb(mp,d,l) + moia(mp,-funcs.rc_2ndstg_p(d,l));
funcs.L1 = @(d,h) d./3.68.*coth(3.68.*h./d);
funcs.m1 = @(mp,d,h) mp.*d./(4.4.*h).*tanh(3.68.*h./d);
funcs.m0 = @(mp,d,h) mp - funcs.m1(mp,d,h);
funcs.l1 = @(d,h) -d./7.36.*csch(7.36.*h./d);
funcs.l0 = @(mp,d,h) (mp.*(h/2-d.^2/(8.*h))-(funcs.l1(d,h)+funcs.L1(d,h)).*funcs.m1(mp,d,h))./funcs.m0(mp,d,h);
% first stage liquid propellant - end burner solid cylinder
funcs.rc_1ststg_p = @(d,l) [-l/2; 0; 0];
funcs.moic_1ststg_p_rgb = @(mp,d,l) mp*diag([1/8*d^2 1/16*d^2+1/12*l^2 1/16*d^2+1/12*l^2]);
funcs.moit_1ststg_p_rgb = @(mp,d,l) funcs.moic_2ndstg_p_rgb(mp,d,l) + moia(mp,-funcs.rc_2ndstg_p(d,l));
end
