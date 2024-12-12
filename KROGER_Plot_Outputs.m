% conditions.figures_flag = 'On';
% % conditions.figures_flag = 'Off';
%
% if strcmp(conditions.figures_flag,'On')


% Clear all the figures to be used
for i=1:7
    figure(i)
    clf
    hold on
end
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%% Plot Formation enthalpies vs EF at 300 K
figure(1)
EF = 0:conditions.EgRef/100:conditions.EgRef;
dH_EF0 = defects.cs_dHo - defects.cs_dm*conditions.muT_equilibrium(1,:)' ;
for i=1:defects.num_chargestates
    f1line(i) = plot(EF, (dH_EF0(i) + defects.cs_charge(i)*EF),'DisplayName', defects.chargestate_names(i));
    f1line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f1line(i).DisplayName},size(f1line(i).XData)));
end

legend show
title(strcat('300 K Formation Enthalpies of Defects for \mu_{Cd}=',num2str(conditions.muT_equilibrium(1,1)),' and \mu_{Te}=',num2str(conditions.muT_equilibrium(1,2))))
xlabel('EF (eV)')
ylabel('Formation Enthalpy (eV)')
datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover

clear EF dH_EF0 i datacursormode f1line
%%%%%%%%%%% end formation enthalpy plotting



%% plot total defect concentrations - these will be same for equilib and quenched except for n and p
figure(2)
f2line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n equilib','MarkerEdgeColor',"r",'MarkerFaceColor',"r");  % full circles for equilib, open for fullquench
f2line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p equilib','MarkerEdgeColor',"b",'MarkerFaceColor',"b");
if conditions.sth_flag==1
    f2line(3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth1),'o','DisplayName','sth1 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
    f2line(4) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth2),'o','DisplayName','sth2 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
    f2line(5) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    f2line(6) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
    f2line(7) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");
    f2line(8) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth1),'o','DisplayName','sth1 quench','MarkerEdgeColor',"g");
    f2line(9) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth2),'o','DisplayName','sth2 quench','MarkerEdgeColor',"g");

    for i = 1:9
        f2line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i).DisplayName},size(f2line(i).XData)));
    end
    legend('n equilib','p equilib', 'sth1 equilib', 'sth2 equilib', '|Nd-Na|', 'n quench', 'p quench',  'sth1 quench', 'sth2 quench' )


    for i=1:defects.num_defects
        if max(log10(equilib_dark_sol.defects(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f2line(i+9) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
            f2line(i+9).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i+9).DisplayName},size(f2line(i+9).XData)));
        end
    end

elseif conditions.sth_flag==0
    f2line(3) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    f2line(4) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
    f2line(5) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");

    for i = 1:5
        f2line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i).DisplayName},size(f2line(i).XData)));
    end
    legend('n equilib','p equilib', '|Nd-Na|', 'n quench', 'p quench' )

    for i=1:defects.num_defects
        if max(log10(equilib_dark_sol.defects(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f2line(i+5) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
            f2line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i+5).DisplayName},size(f2line(i+5).XData)));
        end
    end
end

legend show
title('Concentrations at Equilibrium Temperature K')
xlabel('Equilibrium T (K)')
ylabel('Log_{10} Concentrations (#/cm3)')
ylim([conditions.plot_log10_min conditions.plot_log10_max])
xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover




%% plot all chargestates for equilibrium
figure(3)
f3line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n equilib','MarkerEdgeColor',"r",'MarkerFaceColor',"r");  % full circles for equilib, open for fullquench
f3line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p equilib','MarkerEdgeColor',"b",'MarkerFaceColor',"b");

if conditions.sth_flag==1
    f3line(3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth1),'o','DisplayName','sth1 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
    f3line(4) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth2),'o','DisplayName','sth2 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
    f3line(5) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    f3line(6) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
    f3line(7) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");
    f3line(8) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth1),'o','DisplayName','sth1 quench','MarkerEdgeColor',"g");
    f3line(9) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth2),'o','DisplayName','sth2 quench','MarkerEdgeColor',"g");

    for i = 1:9
        f3line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i).DisplayName},size(f3line(i).XData)));
    end
    legend('n equilib','p equilib', 'sth1 equilib', 'sth2 equilib', '|Nd-Na|', 'n quench', 'p quench',  'sth1 quench', 'sth2 quench' )


    for i=1:defects.num_chargestates
        if max(log10(equilib_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f3line(i+3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
            f3line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i+3).DisplayName},size(f3line(i+3).XData)));
        end
    end

elseif conditions.sth_flag==0
    f3line(3) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    f3line(4) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
    f3line(5) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");

    for i = 1:5
        f3line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i).DisplayName},size(f3line(i).XData)));
    end
    legend('n equilib','p equilib', '|Nd-Na|', 'n quench', 'p quench')


    for i=1:defects.num_chargestates
        if max(log10(equilib_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f3line(i+5) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
            f3line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i+5).DisplayName},size(f3line(i+5).XData)));
        end
    end





end


legend show
title('Concentrations at Equilibrium Temperature K')
xlabel('Equilibrium T (K)')
ylabel('Log_{10} Concentrations (#/cm3)')
ylim([conditions.plot_log10_min conditions.plot_log10_max])
xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover


%
% % plot all chargestates for quenching - chargestates get changed for quenching
% figure(4)
% f4line(1) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n','MarkerEdgeColor',"r");
% f4line(2) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p','MarkerEdgeColor',"b");
% f4line(3) = plot(conditions.T_equilibrium,log10(abs(fullquench_dark_sol.Nd - fullquench_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
% for i = 1:3
%     f4line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i).DisplayName},size(f4line(i).XData)));
% end
% legend('n','p','|Nd-Na|')
% for i=1:defects.num_chargestates
%     if max(log10(fullquench_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
%         f4line(i+3) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
%         f4line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i+3).DisplayName},size(f4line(i+3).XData)));
%     end
% end
% legend show
% title('Concentrations at Equilibrium Temperature K')
% xlabel('Equilibrium T (K)')
% ylabel('Log_{10} Concentrations (#/cm3)')
% ylim([conditions.plot_log10_min conditions.plot_log10_max])
% xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover



%% Plot the  band diagram vs T Ec, Ev, EFn and EFp vs T
figure(5)

% plot(conditions.T_equilibrium,equilib_dark_sol.EFp,'b-')   % note in dark Efn and Efp will lie over eachother
plot(conditions.T_equilibrium,equilib_dark_sol.EFn,'r-')
% plot(conditions.T_equilibrium,fullquench_dark_sol.EFp,'--b')
plot(conditions.T_equilibrium,fullquench_dark_sol.EFn,'--r')

%% handle the T dependent band edges or constant ones
if strcmp(conditions.T_dep_bands_flag,'On')
    plot(conditions.T_equilibrium,conditions.EcT_equilibrium,'k-')
    plot(conditions.T_equilibrium,conditions.EvT_equilibrium,'k-')
    plot(conditions.T_equilibrium,conditions.EcRef*ones(size(conditions.T_equilibrium)),'k--')
    plot(conditions.T_equilibrium,conditions.EvRef*ones(size(conditions.T_equilibrium)),'k--')
    legend('EF(T_{Equilib})', 'EF(T_{Quench})','Ec(T_{Equilib})','Ev(T_{Equilib})','Ec(T_{Quench})','Ev(T_{Quench})')
elseif strcmp(conditions.T_dep_bands_flag,'Off')
    plot(conditions.T_equilibrium,conditions.EcRef*ones(size(conditions.T_equilibrium)),'k--')
    plot(conditions.T_equilibrium,conditions.EvRef*ones(size(conditions.T_equilibrium)),'k--')
    legend('EF(T_{Equilib})', 'EF(T_{Quench})','Ec','Ev')
else
    error('conditions.T_dep_bands_flag must be on or off')
end

title('Fermi Levels for Equilibrium and Quenching vs T')
% legend('EFp(T_{Equilib})','EFn(T_{Equilib})', 'EFp(T_{Quench})', 'EFn(T_{Quench})','Ec(T_{Equilib})','Ev(T_{Equilib})','Ec(T_{Quench})','Ev(T_{Quench})')
% legend('EF(T_{Equilib})', 'EF(T_{Quench})','Ec(T_{Equilib})','Ev(T_{Equilib})','Ec(T_{Quench})','Ev(T_{Quench})')
xlabel('Equilibrium T (K)')
ylabel('Energies (eV)')





%% Plot the mu values for the solution
figure(6)

for i=1:defects.numelements
    plot(conditions.T_equilibrium, equilib_dark_sol.mu)
end

f6line(1) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,1),'DisplayName','Cd');
f6line(2) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,2),'DisplayName','Te');
f6line(3) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,3),'DisplayName','N');
f6line(4) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,4),'DisplayName','P');
f6line(5) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,5),'DisplayName','As');
f6line(6) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,6),'DisplayName','Sb');
f6line(7) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,7),'DisplayName','Cu');

for i = 1:7
    f6line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i).DisplayName},size(f6line(i).XData)));
end
title('Chemical Potentials for Elements at Equilibrium and Quenching vs T')
legend('\mu_{Cd}', '\mu_{Te}', '\mu_{N}', '\mu_{P}', '\mu_{As}', '\mu_{Sb}', '\mu_{Cu}');
xlabel('Equilibrium T (K)')
ylabel('Chem Potentials (eV)')
ylim([-12 0])




%% plot the stoichiometry from equilibrium
figure(7)


for i = 1:conditions.num_elements
    f7line(i) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.element_totals(:,i)),'DisplayName',defects.elementnames(i));
    f7line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i).DisplayName},size(f7line(i).XData)));
end


title('Atomic Concentrations')
legend('Cd','Te','N','P','As','Sb','Cu')
xlabel('Equilibrium T (K)')
ylabel('Log_{10} Concentrations (#/cm3)')
ylim([conditions.plot_log10_min conditions.plot_log10_max])





% analyze the distribution of each element across defects

% here just using it for the frozen ones but could extend this for all
% the elements inclduing those with fixed mu

for ii = 1:numel(conditions.indices_of_fixed_elements)

    figure(7+ii)
    clf
    hold on
    [defects_with_element_fraction,~, defects_with_element_names, ~] = contains_element(equilib_dark_sol, defects, conditions.indices_of_fixed_elements(ii));

    for i=1:size(defects_with_element_fraction,2)
        plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i))
        %     eval(strcat("f",int2str(7+ii),"line(i) = plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i),'DisplayName',defects_with_element_names(i));"))
        %     eval(strcat("f",int2str(7+ii),"line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f",int2str(7+ii),"line(i).DisplayName},size(f"int2str(7+ii),"line(i).XData)));"));
    end

    title('Fraction of Element in Defects')
    legend(defects_with_element_names')
    xlabel('Equilibrium T (K)')
    ylabel('Fraction of Element in Each Defect')
    ylim([0 1])

end




%%% clean up memory
clear datacursormode f1line f11line f2line f3line f4line f5 line f6line f7line i
clear num_Ga num_O num_Si num_H num_Fe num_Sn num_Cr num_Ti num_Ir num_Mg num_Ca num_Zn num_Co num_Zr num_Hf num_Ta num_Ge ii defects_with_element_fraction defects_with_element_names
clear all_elements_out all_defects_cell all_chargestates_out all_defects_out all_defects_cell all_defect_headers all_elements_headers f4line all_elements_cell sig_chargestates_cell sig_chargestates_out element_totals_cell element_totals_headers elements_save_fname
clear all_chargestate_headers all_chargestates_cell save_pname sig_chargestate_headers sig_defect_headers
%%%%%%%%%%%%% end Tasks 1 and 2 (equilib and fullquenching)




















