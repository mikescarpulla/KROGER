%% plot 1 : spatial defect conc plot for specific cooling rate

j = 1;   % index number of cooling rate for the spatial defect plot

y_data = zeros(conditions.SQ_num_chardists, defects.num_defects);
n_data = zeros(conditions.SQ_num_chardists, 1);
figure(2)
hold on

for k = 1:defects.num_defects
    % Extract data for the current defect 
    % y_data = SQ_dark_sol.defects(i, j, :, i_new);  % This should be a 1D vector
    for i = 1:conditions.SQ_num_chardists
        y_data(i,k) = SQ_dark_sol.defects(i, j, SQ_dark_sol.l_value_for_defects_cs(i,j), k);
    end
    % y_data = reshape(SQ_dark_sol.defects(:, j, SQ_dark_sol.l_value_for_defects_cs, i_new), [], 1);  % Ensure y_data is a column vector
    % Plot each curve
    f1line(k) = plot(SQ_dark_sol.SQ_chardists, log10(y_data(:,k)),'DisplayName', defects.defect_names(k));
    f1line(k).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f1line(k).DisplayName},size(f1line(k).XData)));
end
for i = 1:conditions.SQ_num_chardists
    n_data(i,j) = SQ_dark_sol.n(i,j,SQ_dark_sol.l_value_for_defects_cs(i,j));
end
f1line = plot(SQ_dark_sol.SQ_chardists, log10(n_data), 'ro', 'DisplayName', 'n');
for i = 1:conditions.SQ_num_chardists
    p_data(i,j) = SQ_dark_sol.p(i,j,SQ_dark_sol.l_value_for_defects_cs(i,j));
end
f1line = plot(SQ_dark_sol.SQ_chardists, log10(p_data), 'bo', 'DisplayName', 'p');

legend show
title('Concentrations of defects vs char distance at T_{rate} =' ,num2str(conditions.SQ_Trates(j)) )
xlabel('Char distance (cm)')
ylabel('Log_{10} Concentrations (#/cm3)')
ylim([conditions.plot_log10_min conditions.plot_log10_max])
% xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% 
% now make plots of the solution.  Will need to customize this some since
% the arrays will be 3D or 4D since they loop over chardists, Trates, and
% temperatures

%% for outputting vectors for plotting in ORIGIN

n_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for n_i = 1:conditions.SQ_num_chardists
    for n_j = 1:conditions.SQ_num_Trates
        n_data_check = SQ_dark_sol.n(:,:,SQ_dark_sol.l_value_for_defects_cs(n_i,n_j));
        n_data_for_contour(n_i,n_j) = n_data_check(n_i,n_j);
    end
end
n_data_for_contour_column_matrix = reshape(n_data_for_contour,[],1);

p_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for p_i = 1:conditions.SQ_num_chardists
    for p_j = 1:conditions.SQ_num_Trates
        p_data_check = SQ_dark_sol.p(:,:,SQ_dark_sol.l_value_for_defects_cs(p_i,p_j));
        p_data_for_contour(p_i,p_j) = p_data_check(p_i,p_j);
    end
end
p_data_for_contour_column_matrix = reshape(p_data_for_contour,[],1);

n_quench_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for nq_i = 1:conditions.SQ_num_chardists
    for nq_j = 1:conditions.SQ_num_Trates
        n_quench_data_check = SQ_dark_sol.n_quench(:,:,SQ_dark_sol.l_value_for_defects_cs(nq_i,nq_j));
        n_quench_data_for_contour(nq_i,nq_j) = n_quench_data_check(nq_i,nq_j);
    end
end
n_quench_data_for_contour_column_matrix = reshape(n_quench_data_for_contour,[],1);

p_quench_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for pq_i = 1:conditions.SQ_num_chardists
    for pq_j = 1:conditions.SQ_num_Trates
        p_quench_data_check = SQ_dark_sol.p_quench(:,:,SQ_dark_sol.l_value_for_defects_cs(pq_i,pq_j));
        p_quench_data_for_contour(pq_i,pq_j) = p_quench_data_check(pq_i,pq_j);
    end
end
p_quench_data_for_contour_column_matrix = reshape(p_quench_data_for_contour,[],1);

%% T_freeze output vectors for plotting in ORIGIN

T_freeze_V_Cd = SQ_dark_sol.Tfreeze_ij_for_defect_k_holder(:,:,1);
T_freeze_V_Cd_column = reshape(T_freeze_V_Cd,[],1);

T_freeze_Cd_i = SQ_dark_sol.Tfreeze_ij_for_defect_k_holder(:,:,4);
T_freeze_Oi_column = reshape(T_freeze_Cd_i,[],1);

% T_freeze_VGai = SQ_dark_sol.Tfreeze_ij_for_defect_k_holder(:,:,15);
% T_freeze_VGai_column = reshape(T_freeze_VGai,[],1);
% 
% T_freeze_V_Ga = SQ_dark_sol.Tfreeze_ij_for_defect_k_holder(:,:,13);
% T_freeze_V_Ga_column = reshape(T_freeze_V_Ga,[],1);
% 
% T_freeze_V_Ga_Sn_Ga = SQ_dark_sol.Tfreeze_ij_for_defect_k_holder(:,:,155);
% T_freeze_V_Ga_Sn_Ga_column = reshape(T_freeze_V_Ga_Sn_Ga,[],1);


%% Concentration output vectors for plotting in ORIGIN

V_Cd_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for h_i = 1:conditions.SQ_num_chardists
    for h_j = 1:conditions.SQ_num_Trates
        h_data_check = SQ_dark_sol.defects(:,:,SQ_dark_sol.l_value_for_defects_cs(h_i,h_j),1);
        V_Cd_data_for_contour(h_i,h_j) = h_data_check(h_i,h_j);
    end
end
V_Cd_data_for_contour_column_matrix = reshape(V_Cd_data_for_contour,[],1);

Cd_i_defect_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates); % Cd_i(Te)
for l_i = 1:conditions.SQ_num_chardists
    for l_j = 1:conditions.SQ_num_Trates
        l_data_check = SQ_dark_sol.defects(:,:,SQ_dark_sol.l_value_for_defects_cs(l_i,l_j),4);
        Cd_i_defect_data_for_contour(l_i,l_j) = l_data_check(l_i,l_j);
    end
end
Cd_i_data_for_contour_column_matrix = reshape(Cd_i_defect_data_for_contour,[],1);

% V_Ga_Sn_Ga_defect_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
% for l_i = 1:conditions.SQ_num_chardists
%     for l_j = 1:conditions.SQ_num_Trates
%         l_data_check = SQ_dark_sol.defects(:,:,SQ_dark_sol.l_value_for_defects_cs(l_i,l_j),155);
%         V_Ga_Sn_Ga_defect_data_for_contour(l_i,l_j) = l_data_check(l_i,l_j);
%     end
% end
% V_Ga_Sn_Ga_data_for_contour_column_matrix = reshape(V_Ga_Sn_Ga_defect_data_for_contour,[],1);
% 
% Oi_defect_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
% for l_i = 1:conditions.SQ_num_chardists
%     for l_j = 1:conditions.SQ_num_Trates
%         l_data_check = SQ_dark_sol.defects(:,:,SQ_dark_sol.l_value_for_defects_cs(l_i,l_j),6);
%         Oi_defect_data_for_contour(l_i,l_j) = l_data_check(l_i,l_j);
%     end
% end
% Oi_data_for_contour_column_matrix = reshape(Oi_defect_data_for_contour,[],1);
% 
% V_Gai_defect_data_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
% for l_i = 1:conditions.SQ_num_chardists
%     for l_j = 1:conditions.SQ_num_Trates
%         l_data_check = SQ_dark_sol.defects(:,:,SQ_dark_sol.l_value_for_defects_cs(l_i,l_j),15);
%         V_Gai_defect_data_for_contour(l_i,l_j) = l_data_check(l_i,l_j);
%     end
% end
% V_Gai_data_for_contour_column_matrix = reshape(V_Gai_defect_data_for_contour,[],1);

%% x and y axis for contour in ORIGIN

xaxis_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for xi = 1:conditions.SQ_num_chardists
    for yi = 1:conditions.SQ_num_Trates
        xaxis_for_contour(xi,yi) = conditions.SQ_chardists(xi);
    end
end
xaxis_for_contour_column_matrix = reshape(xaxis_for_contour,[],1);

yaxis_for_contour = zeros(conditions.SQ_num_chardists,conditions.SQ_num_Trates);
for xi = 1:conditions.SQ_num_chardists
    for yi = 1:conditions.SQ_num_Trates
        yaxis_for_contour(xi,yi) = conditions.SQ_Trates(yi);
    end
end
yaxis_for_contour_column_matrix = reshape(xaxis_for_contour,[],1);