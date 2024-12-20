function scaledIm = getScaledImage(h, chan)


scaledIm = double(h.stack(:,:,chan));
lims = getFigLims(h.fig, chan, scaledIm);

scaledIm = (scaledIm - lims(1))./(lims(2)-lims(1));

scaledIm(scaledIm>1) = 1;
scaledIm(scaledIm<0) = 0;