function [survival] = analysis( B_cells, number_recycled_b_cells, nb_trial_max, nb_max_B_cells, p_mut, p_recycle, t_cell_selection, conc, p_CDR, final_cycles, nb_Ag)
%

%% GC population
% pop_time = mean(number_recycled_b_cells,1);
% 
% figure(); plot(pop_time);
% %title({['Population of GC b cells over time for 2 Ags with mutations only in the CDR']; ['averaged over ', num2str(n_trial_max), ' trials']; [' with proba recycle = ' num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ', concentration = ' num2str(conc) ' and t cell selection rate = ' num2str(t_cell_selection)]});
% title({['Population of GC b cells for 2 Ags']; ['averaged over ' num2str(n_trial_max) ' trials and conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]});
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6);

figure();
x = linspace(2,200,1000);
curve1 = nb_max_B_cells * (4 * (1-p_mut)^2 * p_recycle * t_cell_selection * conc /(1+conc) ).^(x-2);
plot(x, curve1, ':');

for i = 1:nb_trial_max
    for j = 1:final_cycles(i)
        hold on; plot(number_recycled_b_cells(i,1:j)); 
    end
end
title({['Population of GC b cells for 2 Ags']; [' conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6);

%% GC survival rate
survival = 0;
early = 0;
for i = 1:nb_trial_max
    if final_cycles(i) <= 3
        early = early +1;
    end
    if number_recycled_b_cells(i, final_cycles(i)) ~= 0
        survival = survival + 1;
    end
end
survival = survival / nb_trial_max;
early = early / nb_trial_max;
disp(['ratio of GCs that survive ' num2str(survival)]);
disp(['ratio of GCs that stop early ' num2str(early)]);

%% final cycle
figure(); histogram(final_cycles, 'Normalization', 'probability', 'BinWidth', 10);
title('Number of cycles for each trial', 'Fontweight', 'bold');
set(gca,'FontSize',6);

%% mutations
mutations = zeros(nb_trial_max, size(B_cells,2));
for i = 1:nb_trial_max
    mutations(i,:) = squeeze(B_cells(i,:,nb_Ag+3));
end
figure(); histogram(mutations, 'Normalization', 'probability', 'BinWidth', 1);
title('Number of mutations for each trial', 'Fontweight', 'bold');
set(gca,'FontSize',6);
% 
%% energy and affinity over time
% B_cells size is (nb_trial_max, nb_max_B_cells, nb_Ag + 2);
% Shared_mean_energies = mean(B_cells(:,:,nb_Ag +1),3);
% 
% figure();
% histogram(Shared_mean_energies, 'Normalization', 'probability', 'BinWidth', 10);
% title('Number of cycles for each trial', 'Fontweight', 'bold');
% set(gca,'FontSize',6);
% 
% for i = 1:nb_trial_max
%     for j = 1:final_cycles(i)
%         hold on; plot(Shared_mean_energies(i,1:j)); 
%     end
% end



end
