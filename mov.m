function M = mov(firstframe,lastframe)

%figure('Position',[1 1 700 700])
filename = ['c:\cometmov-' date];
M = makemovie(firstframe,lastframe);
movie2avi(M,filename,'keyframe',15,'fps',15,'compression','None')
movie(M,1);
