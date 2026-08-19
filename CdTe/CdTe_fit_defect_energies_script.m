clear all
clc

%%  Select the .mat file database of defect properties and load it %%%%%%%%
% [conditions.defect_db_fname,conditions.defect_db_pname] = uigetfile('*.mat','Select the defect database file you want to use');
conditions.defect_db_pname = "C:\Users\mscar\Documents\MATLAB\Defects\";
conditions.defect_db_fname = "CdTe_Intuon_natives.mat";
% conditions.defect_db_fname = "CdTe_Intuon_vac_int_only.mat";

load(strcat(conditions.defect_db_pname,conditions.defect_db_fname))
% defect_db_name_root = conditions.defect_db_fname(1:(size(conditions.defect_db_fname,2) - 4)    );  % this takes the filename and strips the . and extension - keep it until savenames established


[output_holder] = CdTe_fit_defect_energies(defects);


disp(output_holder)
figure(1)
plot(output_holder(:,1),output_holder(:,3),'g-')
plot(output_holder(:,1),output_holder(:,4),'g-')


