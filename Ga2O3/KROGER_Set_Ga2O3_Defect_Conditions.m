%% set unidentified shallow doping (fixed charge added to charge balance)
%conditions.Nd = 0;   % imposed shallow doping concentrations
%conditions.Nd = 2e16;   % UID in real samples with Si backgroud
% conditions.Nd = 2e17;   % intentional doping
% conditions.Nd = 0;
% conditions.Na = 0;

% set shallow dopants to zero as default
if ~isfield(conditions,'Nd') || isempty(conditions.Nd)
    conditions.Nd = 0;
end

if ~isfield(conditions,'Na') || isempty(conditions.Na)
    % conditions.Na = 2e17;
    conditions.Na = 0;
end
%%% end shallow doping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
% %% set any interstitials that equilibrate at T_quench
% if strcmp(conditions.interstitial_equilibrate_at_Tquench_flag, 'On')
%     cs_is_simple_interstitial_flag = ~any(defects.cs_num_each_site(:,1:5)) && sum(defects.cs_dm,2)==1;  % if the cs uses none of the {Ga1 Ga2 O1 O2 O3} sites, thus uses only interstital sites AND it adds only one atom, it's a simple interstitial.
%     % also add other defects or chargstates explicitly by manually calling
%     % out their index number here, or a different logic condition
%
% elseif strcmp(conditions.interstitial_equilibrate_at_Tquench_flag, 'Off')
%     cs_is_simple_interstitial_flag = zeros(defects.num_chargestates,1);
% else
%     error('Something weird about some chargestates in terms of whether or not they are simple interstitials')
% end
%
%
% for kk = 1:defects.num_chargestates
%     if  cs_is_simple_interstitial_flag(kk)
%         conditions.set_to_zero_for_paraequilibrium_flag = 1;
%     % elseif some other condition based on diffusivity or, it includes small element like Li or H ...
%     end
% end



%% Set any defects with fixed concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this overrides any conditions set for the elements.
% first initialize a COLUMNvector with the right number of entries for the
% current defect variable.  Then you can go in and st individual ones
conditions.fixed_defects = zeros(defects.num_defects,1);
conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);


% one might want to cancel out some defects from the database to not
% include them in a model run
% conditions.fixed_defects_concentrations(6:7) = 0;   % antisites in Intuon CdTe database
% conditions.fixed_defects(6:7)= 1;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% T independent case - each defect has single value for all temperatures

%
% conditions.fixed_defects(169)= 1;
% % conditions.fixed_defects_concentrations(169) = 8.9e16;   % Si_GaI
% %
% %  %% Sn doped Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 3.6e17;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 4.3e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 2.5e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 6.3e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 1.4e16;   % Mg_GaII
% %
% if index==1 || index==5
%     %% Sn doped Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 3.6e17;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 4.3e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 2.5e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 6.3e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 1.4e16;   % Mg_GaII
%
% elseif index==2 || index==6
% %
%     %% Mg doped
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 2e18;   % Mg_GaII
% %
% elseif index==3  || index==7
%     %% Fe doped
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 5e15;   % Mg_GaII

% elseif index==4  || index==8
%
%     %% UID Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 5e15;   % Mg_GaII
%
% end
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T dependent fixed defect values - use this for example after floating the
% element and finding that only 2 defects dominate but their numbers change
% vs T.  The format will be T down the rows and concentrations for each
% defect across.  So to pick out the column vector giving concentration of
% defect #32 at all the temperatures you'd use conditions.fixed_defects_concentrations(:,32)
%
if strcmp(conditions.T_dep_fixed_defect_flag,'On') && size(conditions.fixed_defects_concentrations,2)==1   % if the T dependent flag is on and the list of fixed defects so far is just a column vector with single values
    conditions.fixed_defects_concentrations = ones(size(conditions.T_equilibrium,2),1) * conditions.fixed_defects_concentrations';  % expand out the list so each T can have a distinct value.

    % now we can go in and assign different values for different defects, for
    % example based on values obtained from a prior simuation.  Yes for now specify the fixed defects
    % inside the if loop checking if the flag is on

    %     conditions.fixed_defects(15)= ones(1,1);  % VGa,ic in the 3/21/24 database
    %% this one for 2100:-20:500
    % conditions.fixed_defects_concentrations(:,15) = [2.26508E+11; 2.01586E+11; 1.85099E+11; 1.76221E+11; 1.74721E+11; 1.80022E+11; 1.96051E+11; 2.218E+11; 2.61093E+11; 3.18136E+11; 3.96893E+11; 5.12249E+11; 6.7033E+11; 8.90973E+11; 1.1998E+12; 1.634E+12; 2.24789E+12; 3.12142E+12; 4.37307E+12; 6.17981E+12; 8.80822E+12; 1.26631E+13; 1.83646E+13; 2.68714E+13; 3.96789E+13; 5.91426E+13; 8.90082E+13; 1.35292E+14; 2.07752E+14; 3.22372E+14; 5.05573E+14; 8.01379E+14; 1.28E+15; 2.08E+15; 3.38E+15; 5.54E+15; 9.03E+15; 1.45E+16; 2.27E+16; 3.37E+16; 4.66E+16; 6.12E+16; 7.31E+16; 8.27E+16; 8.85E+16; 9.05E+16; 8.88E+16; 8.41E+16; 7.71E+16; 6.87E+16; 5.98E+16; 5.03E+16; 4.14E+16; 3.31E+16; 2.60E+16; 1.98E+16; 1.46E+16; 1.07E+16; 7.60E+15; 5.11E+15; 3.43E+15; 2.20E+15; 1.36E+15; 8.11042E+14; 4.64816E+14; 2.55007E+14; 1.34026E+14; 6.59593E+13; 3.22555E+13; 1.47954E+13; 6.45521E+12; 2.69887E+12; 8.25486E+11; 3.0156E+11; 1.10664E+11; 34393838843; 10298916679; 2931857896; 768916197.3; 183911746.7; 39640811.58];
    %% this one for 1500:-100:500
    %     conditions.fixed_defects_concentrations(:,15) = [5.05573E+14; 5.54E+15; 4.66E+16; 9.05E+16; 5.98E+16; 1.98E+16; 3.43E+15; 2.55007E+14; 6.45521E+12];

end   %end the if loop checking if T dependent fixed defects are enabled











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% logic checks on fixed defetcs and sertting default to none fixed
% default is none, catch this.  Note that if elseif loops exit once one
% condition is satisfied

if ~isfield(conditions,'fixed_defects') || sum(isempty(conditions.fixed_defects))
    conditions.fixed_defects = zeros(defects.num_defects,1);
end

if strcmp(conditions.T_dep_fixed_defect_flag,'Off')
    conditions.fixed_defects_concentrations = conditions.fixed_defects_concentrations.*conditions.fixed_defects;          % if a defect is not fixed but by accident it has a target conc set, just set its target concentration to zero also, using the fixed defects vector as a mask

    if ~isfield(conditions,'fixed_defects_concentrations') || sum(isempty(conditions.fixed_defects_concentrations))>0
        conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);
    end

elseif strcmp(conditions.T_dep_fixed_defect_flag,'On')

end

% applies to any fixed defect scenario
% find number of fixed defects and their indices
conditions.indices_of_fixed_defects = find(conditions.fixed_defects);
conditions.num_fixed_defects = numel(conditions.indices_of_fixed_defects);

if conditions.num_fixed_defects~=0
    for i=1:conditions.num_fixed_defects
        msg = strcat('WARNING: Fixed concentration used for...',defects.defect_names(conditions.indices_of_fixed_defects(i)),'...Be sure you meant to do this.');
        disp(msg)
    end
elseif conditions.num_fixed_defects==0
    disp('No Defects fixed.  Proceeding')
end

%%%  end section on fixed defects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%