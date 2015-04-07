function SM_t = getSysMat_multi_src_pos (pos_source, pos_probe_elem, gridx, gridz)
% pos_probe_elem: [x,z]* # of elements, N*2
% pos_source: [x,z]* # of sources, M*2
% gridx: position of x-planes, column vector
% gridz: position of z-planes, column vector
SM_t = getSysMat (pos_source(1,:), pos_probe_elem, gridx, gridz);
for n_src = 2 : size(pos_source, 1)
    SM = getSysMat (pos_source(n_src,:), pos_probe_elem, gridx, gridz);
    SM_t = [SM_t;SM];
end
    
