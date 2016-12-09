function [muts] = analysis( success, B_cells, number_recycled_b_cells, nb_trial_max, nb_max_B_cells, p_mut, p_recycle, t_cell_selection, conc, p_CDR, final_cycles, nb_Ag, nb_cycle_max)
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
title({['Population of GC b cells for 1 Ag']; [' conc = ' num2str(conc) ' proba CDR = ' num2str(p_CDR)]}, 'Fontweight', 'bold');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',10);

%% GC success rate

success = success / nb_trial_max;
disp(['ratio of GCs that are successful ' num2str(success)]);
% survival = 0;
% early = 0;
% for i = 1:nb_trial_max
%     if final_cycles(i) <= 3
%         early = early +1;
%     end
%     if number_recycled_b_cells(i, final_cycles(i)) ~= 0
%         survival = survival + 1;
%     end
% end
% survival = survival / nb_trial_max;
% early = early / nb_trial_max;
% disp(['ratio of GCs that survive ' num2str(survival)]);
% disp(['ratio of GCs that stop early ' num2str(early)]);

%% final cycle
figure(); histogram(final_cycles, 'Normalization', 'probability', 'BinWidth', 5);
title('Number of cycles until GC reaction stops ', 'Fontweight', 'bold');
set(gca,'FontSize',10);
% final_cycles_talied = zeros(nb_cycle_max,1);
% for i = 1:nb_trial_max
%     final_cycles_talied(final_cycles(i),1) = final_cycles_talied(final_cycles(i),1) +1;
% end
% figure(); bar(final_cycles_talied, 0.8);
% title('Number of cycles until GC reaction stops ', 'Fontweight', 'bold');
% set(gca,'FontSize',10);
%% mutations
gusmuts = zeros(50,1);
total = 0;
for i = 1:nb_trial_max
    if number_recycled_b_cells(i, final_cycles(i)) > 0
        for j = 1:number_recycled_b_cells(i, final_cycles(i))
            k = B_cells(i,j,nb_Ag+3);
            gusmuts(k+1,1) = gusmuts(k+1,1)+1;
            total = total +1;
        end
    end
end
muts = gusmuts/total;

figure(); bar(muts, 0.8);
title('Number of mutations for the B cells that are in the GC at the end ', 'Fontweight', 'bold');
set(gca,'FontSize',10);



% for i = 1:nb_trial_max
%     if nnz(muts(i,:)) < 1536
%         count_trial = count_trial +1;
%     end
%     for j = 1:size(B_cells,2)
%         if muts(i,j) == 0
%             count = count +1;
%             %disp(['B cell ' num2str(i) ' trial and ' num2str(j) ' index B cell']);
%             %a = squeeze(B_cells(i,j,:))
%         end
%     end
% end
% disp(['count ' num2str(count)]);
% disp(['count trial ' num2str(count_trial)]);
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
