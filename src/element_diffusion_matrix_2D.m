function [K]=element_diffusion_matrix_2D(v,dof_el,gauss,J,L_el)

% Element diffusion matrix
n_gauss=length(gauss);
K=zeros(dof_el,dof_el);
for i=1:dof_el
    for j=1:dof_el
        for n=1:n_gauss
            K(i,j)=K(i,j)+v*((gauss(n).dW(i).dW'.*2./(L_el))*(gauss(n).dN(j).dN.*2./(L_el')))*gauss(n).w*J;
        end
    end
end

end