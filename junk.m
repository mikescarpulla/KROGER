Emig_copy = defects.cs_Emig;

for d = 1:defects.num_defects
    cs_index = defects.cs_indices_lo(d):defects.cs_indices_hi(d);
    num_cs = numel(cs_index);
    defects.cs_Emig(cs_index) = Emig_copy(d)*ones(num_cs,1);
end