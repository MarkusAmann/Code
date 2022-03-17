function [sol_mesh_t1, sol_mesh_t2, split_idx] = transform_t1(sol_mesh,delta_t1_opt)
    tf = sol_mesh(end);
    split_idx = [find(diff(sol_mesh)==0) find(diff(sol_mesh)==0)+1];
    sol_mesh_t1 = delta_t1_opt*sol_mesh(1:split_idx(1));
    sol_mesh_t2 = sol_mesh(split_idx(2):end);
    sol_mesh_t2_zero_shift = sol_mesh_t2 - sol_mesh_t2(1);
    delta_2 = (tf-sol_mesh_t1(end))/sol_mesh_t2_zero_shift(end);
    sol_mesh_t2_zero_shift_transformed = delta_2*sol_mesh_t2_zero_shift;
    sol_mesh_t2 = sol_mesh_t2_zero_shift_transformed + sol_mesh_t1(end);
end