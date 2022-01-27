//=============================================================================
//
//                 MISSION ANALYSIS : Project TOLOSAT 
//
//=============================================================================

xdel(winsid())
clear
CL_init()


// =====================================================
// SELECT LYDANE PROPAGATION EPOCHS
// =====================================================
// Initial mean keplerian parameters corresponding to the Polar orbit
disp("Initiating mean keplerian parameters");
sma = CL_dataGet("body.Earth.eqRad")+500e3      // Semi-major axis [m]
ecc = 2.e-3                                     // Eccentricity [-]
inc = CL_op_ssoJ2("i", sma, ecc)                // Inclination [rad]
pom = %pi/2;                                    // Argument of perigee [rad]
mlh = 6;                                        // MLTAN [hours]
cjd0 = CL_dat_cal2cjd(2024,01,02,12,0,0);       // Julian date initial time [Julian days]
gom = CL_op_locTime(cjd0, "mlh", mlh, "ra");    // Right ascension (longitude) of the ascending node (RAAN) [rad]
anm = 0;                                        // Mean anomaly [rad]
lydane_kep_ini = [sma;ecc;inc;pom;gom;anm]

lydane_duration=5;
lydane_timestep=30/86400;
lydane_cjd=cjd0+(0:lydane_timestep:lydane_duration);


// =====================================================
// ANALYTICAL PROPAGATION WITH LYDANE
// =====================================================

lydane_kep_results = CL_ex_propagate("lydlp","kep",cjd0,lydane_kep_ini,lydane_cjd,'m');


// Propagate the orbit

// The lydane model is the most complete analytical model of celestlab
// Here, mean parameters are chosen over osculating parameters as they have a better physical significance 
// https://space.stackexchange.com/questions/14731/nuances-of-the-terms-mean-osculating-keplerian-orbital-elements


//// Compute the true and eccentric anomalies from M and ecc
//E_ana = CL_kp_M2E(ecc_ana,M_ana);
//v_ana = CL_kp_M2v(ecc_ana,M_ana);
//
//// Compute the MLTAN fron the RAAN and date
//mltan_ana = CL_op_locTime(cjd, "ra", RAAN_ana, "mlh");
//
//// Compute the evolution of the orbital period
//mm = CL_kp_params('mm',sma_ana);                // Mean motion [rad]
//per = CL_kp_params('per',sma_ana);              // Orbital period [s]
//// /!\ Mean parameter are needed to accurately compute the orbital period.

lydane_sma = lydane_kep_results(1,:);                         // Semi-major axis [m]
lydane_ecc = lydane_kep_results(2,:);                         // Eccentricity [-]
lydane_inc = lydane_kep_results(3,:);                         // Inclination [rad]
lydane_pom = lydane_kep_results(4,:);                         // Argument of perigee [rad]
lydane_RAAN = lydane_kep_results(5,:);                        // Right ascension of the ascending node (RAAN) [rad]
lydane_M = lydane_kep_results(6,:);                           // Mean anomaly [rad]

[sat_pos_tmp,sat_vel_tmp] = CL_oe_kep2car(lydane_kep_results);

// ------------
// HYPOTHESES
// ------------
// Date/time of frame supposed in TREF time scale: 
t0 = CL_dat_cal2cjd(2012,12,1,6,0,0);

// Longitude defining the X axis: 
lon0 = -15 * %pi/180; 

// Orbital elements (sma, ecc, inc, argp, raan, mean anomaly): 
kep_launch = [7000.e3; 0.1; 1; 0; 0.1; 0];

// Note that all frames are fixed with respect to each other
// (All velocities are relative to ECI)  

// Conversion to position and velocity:
[pos_launch, vel_launch] = CL_oe_kep2car(kep_launch); 

// Frame transformation matrix: ECF to "launch frame" at t0:
M1 = CL_rot_angles2matrix(3, lon0); 

// Frame transformation matrix: "ECF" to "ECI" at t0:
M2 = CL_fr_convertMat("ECF", "ECI", t0); 

// Composition of frame transformations (no relative angular velocities): 
// "launch frame" -> ECF followed by:  ECF -> ECI
[M, omega] = CL_rot_compose(M1, [0;0;0], -1, M2, [0;0;0], 1); 

// omega is [0;0;0], but we can still use: 
[pos_eci, vel_eci] = CL_rot_pvConvert(pos_launch, vel_launch, M, omega); 

// Or:
// M = M2*M1'
// pos_eci = M * pos_launch 
// vel_eci = M * vel_launch 

// Convert to orbital elements: 
kep_eci = CL_oe_car2kep(pos_eci, vel_eci); 

disp(kep_eci);



// =====================================================
// PLOTS OF THE ORBIT EVOLUTION
// =====================================================

scf(1);
CL_plot_earthMap();


scf(2);
subplot(231)
plot(lydane_cjd-cjd0,lydane_sma- %CL_eqRad)
title('Altitude')
xlabel('Elapsed days since launch')
ylabel('Altitude (m)')
CL_g_stdaxes();

subplot(232)
plot(lydane_cjd-cjd0,lydane_inc*%CL_rad2deg)
title('Inclination')
xlabel('Elapsed days since launch')
ylabel('Inclination (deg)')
CL_g_stdaxes();

subplot(233)
plot(lydane_cjd-cjd0,lydane_ecc)
title('Eccentricity')
xlabel('Elapsed days since launch')
ylabel('Eccentricity')
CL_g_stdaxes();

subplot(234);
plot(lydane_cjd-cjd0,lydane_pom*%CL_rad2deg)
title('Argument of perigee')
xlabel('Elapsed days since launch')
ylabel('Argument of perigee (deg)')
CL_g_stdaxes();

subplot(235);
plot(lydane_cjd-cjd0,lydane_RAAN*%CL_rad2deg)
title('RAAN')
xlabel('Elapsed days since launch')
ylabel('RAAN (deg)')
CL_g_stdaxes();

subplot(236);
plot(lydane_cjd-cjd0,lydane_M*%CL_rad2deg)
title('Mean Anomaly')
ylabel('Mean Anomaly (deg)')
xlabel('Elapsed days since launch')
CL_g_stdaxes();

scf(2).figure_size=[2000,1000];
deletefile('orbit_evolution_lydane.png');
xs2png(2,'orbit_evolution_lydane.png');


