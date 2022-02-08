function [N] = f_N_plot_2D(csi,eta,p)

% Shape functions for plot
for i=1:length(csi)
    for j=1:length(eta)
        switch p
            
            case 1
                N(1).N(i,j)=1/4*(1-csi(i)).*(1-eta(j));
                N(2).N(i,j)=1/4*(1+csi(i)).*(1-eta(j));
                N(3).N(i,j)=1/4*(1+csi(i)).*(1+eta(j));
                N(4).N(i,j)=1/4*(1-csi(i)).*(1+eta(j));
                
            case 2
                N(1).N(i,j)=+1/4*(1-csi(i))*(1-eta(j))*csi(i)*eta(j);
                N(2).N(i,j)=-1/4*(1+csi(i))*(1-eta(j))*csi(i)*eta(j);
                N(3).N(i,j)=+1/4*(1+csi(i))*(1+eta(j))*csi(i)*eta(j);
                N(4).N(i,j)=-1/4*(1-csi(i))*(1+eta(j))*csi(i)*eta(j);
                N(5).N(i,j)=-1/2*(1-eta(j))*eta(j)*(1-csi(i)^2);
                N(6).N(i,j)=+1/2*(1+csi(i))*csi(i)*(1-eta(j)^2);
                N(7).N(i,j)=+1/2*(1+eta(j))*eta(j)*(1-csi(i)^2);
                N(8).N(i,j)=-1/2*(1-csi(i))*csi(i)*(1-eta(j)^2);
                N(9).N(i,j)=+1*(1-csi(i)^2)*(1-eta(j)^2);
                
        end
    end
end

end