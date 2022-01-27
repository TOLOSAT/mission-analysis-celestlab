load("big_savefile.dat");

duration = t($)-t(1);
[mean_passes_per_day, avg_duration_s, all_passes] = IridiumPassStatistics(ordered_pass_list, duration);

printf(' Statistics : \n    Mean Pass Duration (s) : %f \n', avg_duration_s);
printf('    Mean Passes Per Day : %f \n', mean_passes_per_day);
printf('    Mean Communication Time Per Day (min) : %f \n', mean_passes_per_day*avg_duration_s/60);

// RAANs histogram data
passRAANs = zeros(size(all_passes));
for p=1:length(all_passes)
    passRAANs(p) = raans(t==all_passes(p).start_date_cjd);
end

scf()
histplot(20, passRAANs*%CL_rad2deg, normalization=%f)
xlabel('RAAN (deg)')
ylabel('Number of passes')
title('Pass rate dependency on RAAN')
CL_g_stdaxes()
