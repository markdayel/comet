function M = makemovie(firstframe,lastframe)

frames = (lastframe-firstframe)+1;

M=moviein(frames);
set(gca, 'NextPlot','replacechildren')

for j=firstframe:lastframe
    [filename, err] = sprintf ( 'nodes%05i.txt', j );
    plotnodes(readnodes(filename));
    title(['Frame ',int2str(j),' of ',int2str(lastframe)],'Color','W');
    M(:,((j-firstframe)+1)) = getframe;
end