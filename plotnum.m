function p = plotnum(num)

[filename, err] = sprintf ( 'nodes%05i.txt', num );
plotnodes(readnodes(filename));