function [ratio_stoich, element_totals] = calc_stoich(solution, conditions, defects)

if strcmp(conditions.stoich_flag,'CdTe')
    [ratio_stoich, element_totals] = CdTe_stoich(solution, conditions, defects);

elseif strcmp(conditions.stoich_flag,'Ga2O3')
    [ratio_stoich, element_totals] = Ga2O3_stoich(solution, conditions, defects);

elseif strcmp(conditions.stoich_flag,'GaN')
    [ratio_stoich, element_totals] = GaN_stoich(solution, conditions, defects);
else
    ('Must create and make an eleif option here for the X_stoich() function that computes stoichiometry for a new material');
end

end

