defects.chargestatenames = defects.chargestatenames(indices);
defects.ID = defects.ID(indices);
defects.charge = defects.charge(indices);
defects.dHo = defects.dHo(indices);
defects.dm = defects.dm(indices,:);
defects.sites = defects.sites(indices,:);
defects.pre = defects.pre(indices);
defects.gq_Tfreeze_chargestates = defects.gq_Tfreeze_chargestates(indices);
defects.gq_isfrozen_allT_chargestates = defects.gq_isfrozen_allT_chargestates(indices);
defects.gq_isNOTfrozen_allT_chargestates_index = defects.gq_isNOTfrozen_allT_chargestates_index(indices);

index = [5 6 12 13];
defects.defectnames = defects.defectnames(index);
defects.numchargestates_per_defect = defects.numchargestates_per_defect(index);
defects.gq_Tfreeze_defects = defects.gq_Tfreeze_defects(index);
defects.cs_indices_hi = defects.cs_indices_hi(index);
defects.cs_indices_lo = defects.cs_indices_lo(index);
defects.gq_isfrozen_allT_defects = defects.gq_isfrozen_allT_defects(index);
defects.gq_isfrozen_allT_concdefects = defects.gq_isfrozen_allT_concdefects(index);
defects.gq_isNOTfrozen_allT_defects_index = defects.gq_isNOTfrozen_allT_defects_index(index);

1 6 9 12
