function []=plot_geometry_2D(x,y,x_t,y_t,L_x,L_y,n_gauss,L_el_x,L_el_y,x_gauss,y_gauss,p)

% Plot of geometry

n_el_x=(length(x)-1)/p;
n_el_y=(length(y)-1)/p;

% Screen dimensions
scrsz=get(0,'ScreenSize');
bar=64;

figure('Color',[1 1 1],'Position',[0 0 scrsz(3) (scrsz(4)-bar)])
axes('FontSize',14)

for n=1:n_el_x+1
    plot([x_t(n),x_t(n)],[y_t(1),y_t(n_el_y+1)],'k','LineWidth',2)
    hold on
end

for n=1:n_el_y+1
    plot([x_t(1),x_t(n_el_x+1)],[y_t(n),y_t(n)],'k','LineWidth',2)
end

for i=1:n_el_y
    for j=1:n_el_x
        text(x_t(1)+L_el_x/2+L_el_x*(j-1)+L_el_x/20,y_t(1)+L_el_y*(n_el_y-1/2)-L_el_y*(i-1)+L_el_y/20,num2str((i-1)*n_el_x+j),'Color','k','FontSize',14)
    end
end

n=1;
for i=1:length(y)
    for j=1:length(x)
        plot(x(j),y(end-i+1),'bo','LineWidth',3)
        text(x(j)+L_el_x/10,y(end-i+1)-L_el_y/10,num2str(n),'Color','b')
        n=n+1;
    end
end

n=1;
for i=1:length(y)
    for j=1:length(x)
        if i==1 || mod(i-1,p)==0
            plot(x(j),y(end-i+1),'go','LineWidth',3,'MarkerSize',15)
            text(x(j)+L_el_x/10,y(end-i+1)+L_el_y/10,num2str(n),'Color','g')
            n=n+1;
        else
            if j<=length(x_t)
                plot(x_t(j),y(end-i+1),'go','LineWidth',3,'MarkerSize',15)
                text(x_t(j)+L_el_x/10,y(end-i+1)+L_el_y/10,num2str(n),'Color','g')
                n=n+1;
            end
        end
    end
end

for i=1:n_el_y
    for j=1:n_el_x
        for n=1:n_gauss
            plot(x_t(1)+L_el_x*(j-1)+x_gauss(n),y_t(1)+L_el_y*(i-1)+y_gauss(n),'ro','LineWidth',3)
        end
    end
end

hold off
title('Geometry','FontSize',14)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
grid off
xlim([x(1)-L_x/10,x(end)+L_x/10])
ylim([y(1)-L_y/10,y(end)+L_y/10])

end