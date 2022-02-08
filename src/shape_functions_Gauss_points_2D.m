function [gauss]=shape_functions_Gauss_points_2D(gauss,p)

% Computation of shape functions (and derivatives) at Gauss points
n_gauss=length(gauss);
for n=1:n_gauss
    gauss(n).N=f_N_2D(gauss(n).csi,gauss(n).eta,p);
    gauss(n).dN=f_dN_2D(gauss(n).csi,gauss(n).eta,p);
end

end