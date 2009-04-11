function output=plotnodes(nodesdata)

colordef black
figure(1);
set(1,'Color',[0,0,0])

s = size(nodesdata,1)

%xyz=nodesdata(:,1:3)
x=nodesdata(:,1);
y=nodesdata(:,2);
z=nodesdata(:,3);

polymer=nodesdata(:,4);

rgb=nodesdata(:,5:7);
%g=nodesdata(:,6);
%b=nodesdata(:,7);

%mycmap = get(1,'Colormap');

%fscatter3(x,y,z,rgb);
scatter3(x,y,z,2,rgb);
%plot3(x,y,z,'.','Color',rgb);

%plot3(x(1),y(1),z(1),'Color',rgb(1,:));

%hold on
%for i=2:s
%       plot3(x(i),y(i),z(i),'Color',rgb(i,:));
%end
%hold off

%plot3(x(1:5:s),y(1:5:s),z(1:5:s),'.r');
%hold on
%plot3(x(2:5:s),y(2:5:s),z(2:5:s),'.g');
%plot3(x(3:5:s),y(3:5:s),z(3:5:s),'.b');
%plot3(x(4:5:s),y(4:5:s),z(4:5:s),'.c');
%plot3(x(5:5:s),y(5:5:s),z(5:5:s),'.m');
%hold off

%axis off
%axis([-10 +10 -10 +10 -10 +10])
%set(gca,'Color',[0 0 0]);
%set(gca,'XColor',[1 1 1]);
%set(gca,'YColor',[1 1 1]);
%set(gca,'ZColor',[1 1 1]);
axis([-3 +3 -3 +3 -3 +3])
xlabel('x');
ylabel('y');
zlabel('z');
axis vis3d
grid off
%axis square
