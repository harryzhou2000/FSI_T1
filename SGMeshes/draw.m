


%%
xs = p(1,:);
ys = p(2,:);
ts = t(1:3,:);


for it =  1:size(ts,2)
    patch(xs(ts(:,it)),ys(ts(:,it)),it); 
    
end


%% 
hold on;
for it = 1:size(e,2)
    plot(xs(e(1:2,it)),ys(e(1:2,it)));
    text(xs(e(1,it)),ys(e(1,it)),num2str(it));
end
hold off;


%%
file = fopen('mesh1.txt','w');
fprintf(file, '%d\t%d\t%d\n',size(p,2),size(ts,2),size(e,2));
fprintf(file,'P\n');
for it =  1:size(xs,2)
    fprintf(file,'%.7e\t%.7e\n',xs(it),ys(it));
end
fprintf(file,'T\n');
for it =  1:size(ts,2)
    fprintf(file,'%d\t%d\t%d\n',ts(1,it),ts(2,it),ts(3,it));
end
fprintf(file,'ED\n');
for it = 1:size(e,2)
    fprintf(file,'%d\t%d\n',e(1,it),e(2,it));
end
fclose('all');


%%
ttest = ts(:,1);
patch(xs(ttest),ys(ttest),1);
for it = 1:3
    text(xs(ttest(it)),ys(ttest(it)),num2str(it));
end
x1 = xs(ttest(1))-xs(ttest(2));
x2 = xs(ttest(2))-xs(ttest(3));
y1 = ys(ttest(1))-ys(ttest(2));
y2 = ys(ttest(2))-ys(ttest(3));
area = x2*y1-x1*y2;


%%



