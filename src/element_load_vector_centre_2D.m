function [f]=element_load_vector_centre_2D(s_c,dof_el_v,gauss,J)

% Element load vector
n_gauss=length(gauss);
f=zeros(dof_el_v,1);
for i=1:dof_el_v
    for n=1:n_gauss
        f(i,1)=f(i,1)+(s_c*gauss(n).W(i).W)*gauss(n).w*J;
    end
end

end