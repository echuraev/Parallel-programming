clear h;

h.dir = './res/';
h.n = 1999;
h.outfile = 'swe.gif';
h.drawstep = 5;

h.xstep = 1;
h.ystep = 1;
h.xmin = 0;
h.ymin = 0;
h.xmax = 120 - h.xstep;
h.ymax = 120 - h.ystep;

h.x = [h.xmin:h.xstep:h.xmax];
h.y = [h.ymin:h.ystep:h.ymax];
h.arr = importdata([h.dir num2str(0) '.csv']);

h.fig = figure('Position', [0,500,1200,900], 'PaperPositionMode', 'auto');
h.f = surf(h.x, h.y, h.arr);
colormap(winter);
axis([h.xmin h.xmax h.ymin h.ymax -1 3])
for i = 1:h.drawstep:h.n
    h.arr = importdata([h.dir num2str(i) '.csv']);
    set(h.f, 'zdata', h.arr);
    drawnow;
    
    h.frame = getframe(h.fig);
    h.im = frame2im(h.frame);
    [h.gif, h.map] = rgb2ind(h.im, 256);
    if (i == 1)
        imwrite(h.gif, h.map, h.outfile, 'DelayTime', 0.3, 'LoopCount', inf);
    else
        imwrite(h.gif, h.map, h.outfile, 'WriteMode', 'append');
    end;
end

clear h;