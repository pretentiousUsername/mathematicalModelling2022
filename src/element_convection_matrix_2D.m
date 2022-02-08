function [C]=element_convection_matrix_2D(a,dof_el,gauss,J,L_el)

% Element convection matrix
n_gauss=length(gauss);
C=zeros(dof_el,dof_el);
for i=1:dof_el
    for j=1:dof_el
        for n=1:n_gauss
            C(i,j)=C(i,j)+(gauss(n).W(i).W*(a*(gauss(n).dN(j).dN.*2./L_el')))*gauss(n).w*J;
        end
    end
end

end