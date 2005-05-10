function output = tracknodes(startfilenum, endfilenum)

startnode = 1;
endnode = 200;

for j=startfilenum:endfilenum
    [filename, err] = sprintf ( 'nodes%05i.txt', j );
    nodes=readnodes(filename);    
    indices = find( startnode <= nodes(:,1) & nodes(:,1) <= endnode);
    output(:,:,1+(j-startfilenum)) = nodes(indices,2:4);
end

cmap = colormap(hsv((endnode-startnode)));

hold on

for j=1:(endnode-startnode)
        x = output(j,1,:);
    y = output(j,2,:);
    z = output(j,3,:);
    %scatter3(x,y,z,'.');
    %for i=1:(endfilenum-startfilenum)
        plot3(x(:),y(:),z(:),'-','Color',cmap(j,:));
    %end
end

axis([-4 +4 -4 +4 -4 +4])
%axis square
axis vis3d
xlabel('x');
ylabel('y');
zlabel('z');
whitebg('black');
grid off
%plot3(output(:,1,:),output(:,2,:),output(:,3,:));

axis off
set(gcf,'Color','black');
 set(gcf,'Units','pixels','Position',[50 50 651 652])


M=moviein(90);
camtarget([0 0 0]);
set(gca, 'NextPlot','replacechildren')
for j=1:90
view(j*4,10)
M(:,j) = getframe;
end

movie(M,4);

%filename = ['c:\nodetracking9-' date];
%movie2avi(M,filename,'keyframe',15,'fps',15,'compression','None')

