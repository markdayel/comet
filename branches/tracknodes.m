function output = tracknodes(startfilenum, endfilenum)

startnode = 1;
endnode = 200;


[filename, err] = sprintf ( 'nodes%05i.txt', 1 );
nodes=readnodes(filename); 

%get the nucleator shape
indices = find( 0 == nodes(:,1));
nucleator = nodes(indices,2:4);
    
for j=startfilenum:endfilenum
    [filename, err] = sprintf ( 'nodes%05i.txt', j );
    nodes=readnodes(filename);    
    indices = find( startnode <= nodes(:,1) & nodes(:,1) <= endnode);
    output(:,:,1+(j-startfilenum)) = nodes(indices,2:4);
end

cmap = colormap(hsv((endnode-startnode)));

hold on

plot3(nucleator(:,1),nucleator(:,2),nucleator(:,3),'.w','MarkerSize',1);

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
 set(gcf,'Units','pixels','Position',[50 50 452 550])


M=moviein(640);
camtarget([0 0 0]);
set(gca,'projection','perspective')
set(gca, 'NextPlot','replacechildren')
framecount = 1;
for j=0:0.025:(2*3.141)
    set(gca,'cameraposition',[20*cos(j),20*sin(j),0])
    set(gca,'cameratarget',[0,0,0])
    set(gca,'cameraupvector',[0 0 1])
    set(gca,'cameraviewangle',15)
    set(gca,'xlim',[-4 4],'ylim',[-4 4],'zlim',[-4 4])
    box on
    %view(j*4,10)
M(:,framecount) = getframe(gcf);
framecount = framecount + 1;
end



movie(M,2);

filename = ['c:\nodetracking-' date];
movie2avi(M,filename,'keyframe',15,'fps',15,'compression','None')

