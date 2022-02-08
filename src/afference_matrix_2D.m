function [A]=afference_matrix_2D(n_np_x,n_np_y,dof_el,p)

% Afference matrix
switch p
    case 1
        n_el=(n_np_x-1)*(n_np_y-1);
        A=zeros(n_el,dof_el);
        for i=1:n_el
            [r,c]=row_column(i,n_np_x-1);
            A(i,1)=r*n_np_x+1+(c-1);
            A(i,2)=r*n_np_x+2+(c-1);
            A(i,3)=(r-1)*n_np_x+2+(c-1);
            A(i,4)=(r-1)*n_np_x+1+(c-1);
        end
        
    case 2
        n_el_x=(n_np_x-1)/2;
        n_el=(n_np_x-1)/2*(n_np_y-1)/2;
        A=zeros(n_el,dof_el);
        for i=1:n_el
            [r,c]=row_column(i,n_el_x);
            A(i,1)=      r*2*n_np_x+1+2*(c-1);
            A(i,2)=      r*2*n_np_x+3+2*(c-1);
            A(i,3)=  (r-1)*2*n_np_x+3+2*(c-1);
            A(i,4)=  (r-1)*2*n_np_x+1+2*(c-1);
            A(i,5)=      r*2*n_np_x+2+2*(c-1);
            A(i,6)=(r-1/2)*2*n_np_x+3+2*(c-1);
            A(i,7)=  (r-1)*2*n_np_x+2+2*(c-1);
            A(i,8)=(r-1/2)*2*n_np_x+1+2*(c-1);
            A(i,9)=(r-1/2)*2*n_np_x+2+2*(c-1);
        end
end

end