function out = z_loc (I,z_slice)

z=1:z_slice;

opts_e = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Lower', [0, 1, 0], ...    % Lower bounds: a >= 0, b is unconstrained
    'Upper', [65535, z_slice, 65535]);



pd_e=fit(z', I,'gauss1',opts_e);
par_e=coeffvalues(pd_e);

out = par_e(1,2);