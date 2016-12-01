function [ pop_time, survival, total_exit_cells, neutralized, breadth, average_energy, affinity ] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, n_trial_max, e_act, n_cycle_max, p_mut, p_recycle, t_cell_selection, conc, p_CDR)
%

%% GC population
pop_time = mean(number_recycled_b_cells,1);

figure(); plot(pop_time);
%title({['Population of GC b cells over time for 2 Ags with mutations only in the CDR']; ['averaged over ', num2str(n_trial_max), ' trials']; [' with proba recycle = ' num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ', concentration = ' num2str(conc) ' and t cell selection rate = ' num2str(t_cell_selection)]});
title({['Population of GC b cells for 2 Ags']; ['averaged over ' num2str(n_trial_max) ' trials and conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]});
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6);

figure(); 
for i = 1:n_trial_max
    plot(number_recycled_b_cells(i,:)); hold on;
end
title({['Population of GC b cells for 2 Ags']; [' conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6);

%% GC survival rate
survival = 0;
for i = 1:n_trial_max
    if number_recycled_b_cells(i,n_cycle_max -1) ~= 0
        survival = survival + 1;
    end
end
survival = survival / n_trial_max;
disp(['ratio of GCs that survive ' num2str(survival)]);

%% breadth 
% fraction of exit_cells with affinity for both Ags above a threshold for different thresholds
% B_cells = zeros(nb_trial_max, nb_max_B_cells, nb_Ag + 2);
% exit_cells = zeros(n_trial_max, n_cycle_max, n_Ag, floor(n_max_Bcells/4));
% 
% total_exit_cells = mean(number_exit_cells,1);
% thresholds = linspace(e_act-3, e_act+2,5);
% neutralized = zeros(n_trial_max, n_cycle_max, length(thresholds));
% breadth = zeros(n_cycle_max,length(thresholds));
% 
% for t = 1:length(thresholds)
%     for i = 1:size(exit_cells,1) %trial
%         for j = 1:size(exit_cells,2) %cycle
%             for k = 1:size(exit_cells,4) %bcell index
%                    if (exit_cells(i,j,1,k) >= thresholds(t) && exit_cells(i,j,2,k) >= thresholds(t))
%                         neutralized(i,j,t) = neutralized(i,j,t) +1;
%                    end
%             end
%         end
%     end
% end
% 
% neutralized_mean = mean(neutralized,1);
% 
% figure();
% plot(total_exit_cells);
% title({['Number exit cells for 2 Ags']; ['averaged over ' num2str(n_trial_max) ' trials and conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold');
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6);
% 
% figure();
% for t = 1:length(thresholds)
%     plot(neutralized_mean(:,:,t)); hold on;
% end
% title({['Neutralized for 2 Ags']; ['averaged over ' num2str(n_trial_max) ' trials and conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold');
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6);
% legendCell = strcat(strtrim(cellstr(num2str(thresholds(:)))));
% legend(legendCell, 'fontsize',6 , 'Position', [0.75,0.65,0.15,0.25]);
% 
% for c = 1:n_cycle_max
%     if total_exit_cells(c) == 0
%         breadth(c,t) = 0;
%     else
%         for t = 1:length(thresholds) 
%             breadth(c,t) = neutralized_mean(1,c,t);
%         end
%         breadth(c,:) = breadth(c,:)/total_exit_cells(c);
%     end
% end
% 
% figure();
% for t = 1:length(thresholds)
%     plot(breadth(:,t)); hold on;
% end
% title({['Breadth for 2 Ags']; ['averaged over ' num2str(n_trial_max) ' trials and conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold')
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6);
% legend(legendCell, 'fontsize',6, 'Position', [0.75,0.65,0.15,0.25]);
% 
% %% energy and affinity over time
% % exit_cells = zeros(n_trial_max, n_cycle_max, n_Ag, floor(n_max_Bcells/4));
% exit_cells_averaged_trials = mean(exit_cells,1);
% average_energy = zeros(1,n_cycle_max);
% for c = 1:n_cycle_max
%     if total_exit_cells(c) == 0
%         average_energy(c) = 0;
%     else
%         for k = 1:size(exit_cells,4) %sum the average energy per cycle
%             average_energy(c) = average_energy(c) + 0.5*(exit_cells_averaged_trials(1,c,1,k)+exit_cells_averaged_trials(1,c,2,k));
%         end
%         average_energy(c) = (average_energy(c)/total_exit_cells(c))-e_act;
%     end
% end
% 
% figure();
% plot(average_energy);
% title('Energy increase');
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6)
% 
% figure();
% affinity = exp(average_energy);
% plot(affinity);
% title('Affinity');
% xlabel('Number of cycles', 'Fontweight', 'bold');
% set(gca,'FontSize',6)

end
