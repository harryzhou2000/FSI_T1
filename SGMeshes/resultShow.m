

getresult;
hold on
for it = 1:30:size(result,1)
    if isnan(result(it,2))
        x = result((it+1):(it+3),1);
        y = result((it+1):(it+3),2);
        z = result((it+1):(it+3),3);
        if 1
            plot3(x,y,z,'.k');
            %%patch(x,y,z,z);
        end
    end
end
hold off

view(30,45);

[xm,ym]=meshgrid(0:0.05:1,0:0.05:1);
hold on;
surf(xm,ym,-0.25*(xm.^2+ym.^2)+0.0,'LineStyle','none');

hold off;