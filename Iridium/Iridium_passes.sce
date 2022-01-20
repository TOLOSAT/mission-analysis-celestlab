// don't forget to run Fonction_iridium_utile.sce before running this ! 
// you will need Celestlab and CelestlabX

clear all;
CL_init();

//= = = = = = = = = = = = = = = = = = = = = = =  Variables  = = = = = = = = = = = = = = = = = = = = = = =  

// with these settings it takes my computer 4 minutes to run the code
// but running over longer durations yields much more statistically interesting results
time_step = 5/86400;   // in days
duration = 1*24/24; // in days
min_pass_duration_s = 15;
min_pass_size = ceil(15/(time_step*86400));

t0 = CL_dat_cal2cjd(2017,01,16,12,0,0); //CL_dat_cal2cjd(2024,07,02,12,0,0); //year, month day;    // initial time
t = t0 + (0:time_step:duration) ;
freq_emission = 1617.e6;
Delta_doppler_max = 345; // +/- 345Hz/s max Doppler rate (time-derivative of Doppler Shift)
doppler_max = 35000; // Hz

//= = = = = = = = = = = = = = = = = = = = = = =  generer la constellation  = = = = = = = = = = = = = = = = = = = = = = =  

[kepConst, sat_angle_lim, cen_angle_lim, freq_emission_const, freq_reception_const] = generate_constellation("IridiumNext");


//= = = = = = = = = = = = = = = = = = = = = = =  generation satellite  = = = = = = = = = = = = = = = = = = = = = = =  

//  ==>  ==>  ==>  ==>  ==>  ==>  ==>  ==>  Polar orbit 

sma = 6748.e3 //6878.e3   //sma 
ecc = 2.e-3 //0.0253508
inc = 51.6*%CL_deg2rad //97.4*%CL_deg2rad
pom = 1.570796327;
mlh = 6; // MLTAN (hours)
gom = 3.3392258;
anm = 0;
kep_mean_ini = [  sma;     // demi grand axe (m)
ecc;     // excentricité
inc;     // inclinaison (rad)
pom;     // argument du périastre, petit omega (rad)
gom;     // longitude du noeud ascendant raan (rad)  
anm] 

T = 2*%pi*sqrt(sma^3/%CL_mu);

// the value 62.9deg below comes from the Iridium satellite half cone opening angle IIRC
min_el_for_vizi_deg = %CL_rad2deg * acos(sin(62.9*%CL_deg2rad)*kepConst(1,1)/sma); // about 23deg

//= = = = = = = = = = = = = = = = = = = = = = =  generation  = = = = = = = = = = = = = = = = = = = = = = =  

kep_sat = CL_ex_secularJ2(t0,kep_mean_ini,t); 
t_since_periapsis = pmodulo((t-t0)*86400,T);
MA_deg = t_since_periapsis*360/T;
aol = pmodulo(kep_sat(6,:) + kep_sat(4,:), 2*%pi); // mean argument of latitude is the interesting parameter

//convertion to cartesian position and velocity:
[pos_sat, vel_sat] = CL_oe_kep2car(kep_sat);
[pos_sat,vel_sat]= CL_fr_convert("ECI", "ECF", t, pos_sat,[vel_sat]);
[i,j]=size(pos_sat);

//= = = = = = = = = = = = = = = = = = = = = = =  simulation   = = = = = = = = = = = = = = = = = = = = = = =  

// initializations
t1=t0;
t1_old = t1;
[n, m] = size(kepConst); // n=6 orbital elements, m=75 satellites
elevations_over_t = zeros(m,length(t)); //deg
doppler_shifts = zeros(m,length(t)); // HZ
doppler_rates = zeros(m,length(t)); // Hz/s

// constellation propagation
kepConstIridium = zeros(n, m, length(t));
for k=1:m
    kepConstIridium(:,k,:) = CL_ex_secularJ2(t0, kepConst(:,k), t);
end
    
// LOOP ON TIME
// positions have been precalculated
// calcs elevation and doppler shift
for l=1:j // (lowercase L), t = t0+dt*l
    // get params of all iridium sats at tl
    kep = kepConstIridium(:,:,l);
    [pos, vel] = CL_oe_kep2car(kep);
    [pos,vel]= CL_fr_convert("ECI", "ECF", t1,pos,[vel]);

    // compute angles between our sat and the Iridium sats
    sat_angle = CL_vectAngle(CL_dMult(ones(1,m),pos_sat(:,l))-pos, -pos);
    cen_angle = CL_vectAngle(pos_sat(:,l), pos);

    // For each Iridium sat
    // calc doppler shift and elevation of iridium sat seen by our sat
    for k=1:m // k is iridium index
        shift = doppler_shift(pos(:,k), vel(:,k), pos_sat(:,l), vel_sat(:,l), freq_emission);
        el = 90-sat_angle(k)*%CL_rad2deg-cen_angle(k)*%CL_rad2deg;
        elevations_over_t(k,l) = el;
        doppler_shifts(k,l) = shift;
    end 
    t1_old = t1;
    t1=t1+time_step;
end
// second loop : just to compute doppler rates a fortiori using centered numerical derivative scheme (where possible)
for k=1:m
    doppler_rates(k,1) = (doppler_shifts(k,2)-doppler_shifts(k,1))/(time_step*86400);
    for l=2:(j-1)
        doppler_rates(k,l) = (doppler_shifts(k,l+1)-doppler_shifts(k,l-1))/(2*time_step*86400);
    end
    doppler_rates(k,j) = (doppler_shifts(k,j)-doppler_shifts(k,j-1))/(time_step*86400);
end

//= = = = = = = = = = = = = = = = = = = = = = =  exploitation  = = = = = = = = = = = = = = = = = = = = = = =

// visibility of each sat at each time step
is_visible = (elevations_over_t >= min_el_for_vizi_deg & ... 
              abs(doppler_shifts) <= doppler_max & ...
              abs(doppler_rates) <= Delta_doppler_max)*1; // 0s and 1s


// statistics on iridium passes
// 1. pass structure definition
example_pass = struct('start_date_cjd', t(100), 'end_date_cjd', t(106), 'sat_number', 46);
// 2. helper functions
function pass_duration = IridiumPassDuration_s(IridiumPass)
    pass_duration = (IridiumPass.end_date_cjd - IridiumPass.start_date_cjd)*86400;
endfunction
// 3. get all passes
ordered_pass_list = list(); // dimension 1 is sat number, dim 2 is time
for k=1:m
    ordered_pass_list($+1) = list();
end
for l_start = 2:j // loop on pass start time
    for k=1:m // loop on sat number
        if (is_visible(k,l_start) & ~is_visible(k,l_start-1)) // rising edge detected
            // now search for falling edge (or end of data)
            maybe_end_l = l_start;
            while ((maybe_end_l <= j) & is_visible(k,maybe_end_l))
                maybe_end_l = maybe_end_l + 1;
            end
            // now maybe_end_l is either j+1 or the true end index of the pass
            l_end = maybe_end_l - 1;
            // l_start is the first index where the sat is visible. 
            // l_end is the last index where the sat is visible.
            // we ignore passes that are shorter than min_pass_duration
            if (l_end-l_start)*time_step*86400 < min_pass_duration_s
                continue
            end
            // keep this pass if another pass is not already happening
            if sum(is_visible(:,l_start)) < 1.5
                pass = struct('start_date_cjd', t(l_start), 'end_date_cjd', t(l_end), 'sat_number', k);
                ordered_pass_list(k)($+1) = pass;
            end
        end
    end
end
// 4. statistics of passes
N_passes = 0;
avg_duration_s = 0;
max_duration_s = 0;
min_duration_s = 99999999;
total_duration_s = 0;
all_passes = list(); // unordered 1D list
for k=1:m
    for p=1:length(ordered_pass_list(k))
        N_passes = N_passes+1;
        pass = ordered_pass_list(k)(p);
        all_passes($+1) = pass;
        dur = IridiumPassDuration_s(pass);
        avg_duration_s = avg_duration_s + dur;
        max_duration_s = max(max_duration_s, dur);
        min_duration_s = min(min_duration_s, dur);
        total_duration_s = total_duration_s + dur;
    end
end
avg_duration_s = avg_duration_s/N_passes;
mean_passes_per_day = N_passes/duration;
printf('\n\n Statistics : \n    Mean Pass Duration (s) : %f \n', avg_duration_s);
printf('    Max Pass Duration (s) : %f \n', max_duration_s);
printf('    Min Pass Duration (s) : %f \n', min_duration_s);
printf('    Mean Passes Per Day : %f \n', mean_passes_per_day);
printf('    Mean Communication Time Per Day (min) : %f \n', total_duration_s/60/duration);

// compute expected number (in statistical sense) of  iridium satellites
// that satisfy visibility + doppler shift + doppler rate constraints
// as a function of AoL (argument of latitude, = omega + TA)
// objective : show that some geographic areas are more favorable to communication
visi_sat_expected_over_aol = zeros(1, length(t((t-t0)*86400<T)));
daol = time_step*86400*2*%pi/T;
for l=1:length(t((t-t0)*86400<T))
    myaol = aol(1,l);
    visi_sat_expected_over_aol(1,l) = sum(is_visible(:,(myaol <= aol & (aol < (myaol+daol))) ));
end
visi_sat_expected_over_aol = visi_sat_expected_over_aol./(duration.*86400./T);

// = = = = = = = = = = = = = = = = = = = = = = =  visualisation  = = = = = = = = = = = = = = = = = = = = = = =

scf()
plot((t-t0)*24*60, is_visible)
xlabel('elapsed time (min)')
ylabel('Visible satellites')
title('Visible satellites ')
CL_g_stdaxes()
set(gca(),'data_bounds',[-10, -0.2; 10+duration*24*60, 1.2])

// pass analytics 
// plots elevation, doppler shift and doppler rate over time for each interesting pass
// TODO: modify this code so it works if 1 sat has several passes
//for k=1:m
//    if (max(is_visible(k,:)) == 1 & duration < 1) then
//        indices_of_pass = elevations_over_t(k,:) >= min_el_for_vizi_deg;
//        scf()
//        subplot(3,1,1)
//        plot((t(indices_of_pass)'-t0)*1440,elevations_over_t(k,indices_of_pass))
//        xlabel('temps ecoulé (min)')
//        ylabel('elevation [deg]')
//        title('iridium ' + string(k))
//        CL_g_stdaxes()
//    
//        subplot(3,1,2)
//        plot((t(indices_of_pass)'-t0)*1440,doppler_shifts(k,indices_of_pass))
//        xlabel('temps ecoulé (min)')
//        ylabel('doppler shift [Hz]')
//        CL_g_stdaxes()
//        
//        subplot(3,1,3)
//        plot((t(indices_of_pass)'-t0)*1440,doppler_rates(k,indices_of_pass))
//        xlabel('temps ecoulé (min)')
//        ylabel('doppler rate [Hz/s]')
//        CL_g_stdaxes()
//    end
//end


duration_min = duration*24*60
total_comm_time_min = sum(is_visible)*time_step*24*60
comm_availability_ratio = total_comm_time_min/duration_min
// TODO : get individual passes length
scf()
plot(MA_deg, is_visible)
xlabel('mean anomaly (deg)')
ylabel('Visible satellites')
title('Visible satellites ')
CL_g_stdaxes()
set(gca(),'data_bounds',[-5, -0.2; 365, 1.2])
//scf()
//plot(t_since_periapsis/60, is_visible)
//xlabel('time since periapsis [min]')
//ylabel('Visible satellites')
//title('Visible satellites ')
//CL_g_stdaxes()
//set(gca(),'data_bounds',[-5, -0.2; 115, 1.2])
scf()
plot(aol*180/%pi, sum(is_visible,1))
xlabel('argument of latitude [deg]')
ylabel('Visible satellites')
title('Visible satellites ')
CL_g_stdaxes()
set(gca(),'data_bounds',[-5, -0.2; 365, 2.2])

scf()
plot(aol((t-t0)*86400<T)*180/%pi, visi_sat_expected_over_aol, 'x')
xlabel('argument of latitude [deg]')
ylabel('Expected Number of Visible Satellites')
title('Expected Visible Satellites, averaged over XX days, dt=XXs ')
CL_g_stdaxes()
set(gca(),'data_bounds',[-5, -0.2; 365, 1.2])
scf()
plot(t-t0, sum(is_visible, 1))
xlabel('mission elapsed time (days)')
ylabel('Visible satellites')
title('Visible satellites ')
CL_g_stdaxes()

raans = kep_sat(5,:);
smas = kep_sat(1,:);
scf()
plot(t-t0, raans)
xlabel('mission elapsed time (days)')
ylabel('RAAN [rad]')
title('Right Ascension of Ascending node, 10 days at 370km 51.6deg')
CL_g_stdaxes()
scf()
plot(t-t0, smas)
xlabel('mission elapsed time (days)')
ylabel('SMA [km]')
title('SMA, 10 days at 370km 51.6deg')
CL_g_stdaxes()
