function SM = getSysMat (pos_source, pos_probe_elem, gridx, gridz)
% Fereshteh Added: 
% Assumptions: No ray is perpendicular to Z-axis, dx and dz are constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pos_probe_elem: [x,z]* # of elements, N*2
% pos_source: [x,z], 1*2
% gridx: position of x-planes, column vector
% gridz: position of z-planes, column vector
zero_tol = 1e-8;
px = pos_source(1);
pz = pos_source(2);
Nx = length(gridx); % number of planes
Nz = length(gridz); % number of planes
dx = abs(gridx(2) - gridx(1));
dz = abs(gridz(2) - gridz(1));
SM = sparse(size(pos_probe_elem,1), (Nx - 1)*(Nz-1));
for ind_ele = 1:size(pos_probe_elem,1) % for each ray
    ele_x = pos_probe_elem(ind_ele,1);
    ele_z = pos_probe_elem(ind_ele,2);
    d_source_ele = sqrt((ele_x - px)^2 + (ele_z - pz)^2);
    alphaz_1 = (gridz(1) - pz)/(ele_z - pz);
    alphaz_Nz = (gridz(end) - pz)/(ele_z - pz);
    if abs(px - ele_x) > zero_tol    
        alphax_1 = (gridx(1) - px)/(ele_x - px);
        alphax_Nx = (gridx(end) - px)/(ele_x - px);
        alpha_min = max([0,min(alphax_1, alphax_Nx),min(alphaz_1, alphaz_Nz)]);
        alpha_max = min([1,max(alphax_1, alphax_Nx),max(alphaz_1, alphaz_Nz)]);
    else
        alpha_min = max(0, min(alphaz_1, alphaz_Nz));
        alpha_max = min(1, max(alphaz_1, alphaz_Nz));
    end
    
    if (ele_x - px) >= 0
        i_min = ceil(Nx - (gridx(end) - alpha_min*(ele_x - px) - px)/dx);
        i_max = floor(1 + (px + alpha_max*(ele_x - px) - gridx(1))/dx);
    else
        i_min = ceil(Nx - (gridx(end) - alpha_max*(ele_x - px) - px)/dx);
        i_max = floor(1 + (px + alpha_min*(ele_x - px) - gridx(1))/dx); 
    end

    if (ele_z - pz) >= 0
        j_min = ceil(Nz - (gridz(end) - alpha_min*(ele_z - pz) - pz)/dz);
        j_max = floor(1 + (pz + alpha_max*(ele_z - pz) - gridz(1))/dz);
    else
        j_min = ceil(Nz - (gridz(end) - alpha_max*(ele_z - pz) - pz)/dz);
        j_max = floor(1 + (pz + alpha_min*(ele_z - pz) - gridz(1))/dz); 
    end

    alpha_z = (gridz - pz) / (ele_z - pz); % column vector
    if abs(px - ele_x) > zero_tol
        alpha_x = (gridx - px) / (ele_x - px); % column vector
        alpha = unique([alpha_min;alpha_max;alpha_x;alpha_z]); % should be (n+1)*1 vector
    else
        alpha = unique([alpha_min;alpha_max;alpha_z]); % should be (n+1)*1 vector
    end
    alpha (alpha < alpha_min | alpha > alpha_max) = [];
    
    n = length(alpha) - 1;
    
    for m = 1:n
        alpha_mid = (alpha(m) + alpha(m+1))/2;
        ind_x = floor(1 + (px + alpha_mid*(ele_x - px) - gridx(1))/dx);
        ind_z = floor(1 + (pz + alpha_mid*(ele_z - pz) - gridz(1))/dz);
        d_pix = d_source_ele * (alpha(m+1) - alpha(m));
        ind_SM_col = (ind_x - 1)*(Nz - 1) + ind_z;
        SM(ind_ele, ind_SM_col) = d_pix;
    end
    
end
        
        
