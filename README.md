# costMatLinearGuidedTracks
A set of cost matrices for tracking axonal transport using uTrack.

uTrack (http://lccb.hms.harvard.edu/software.html) is a powerful particle tracking framework for Matlab developed by Khuloud Jaqaman et al. [see also http://www.nature.com/nmeth/journal/v5/n8/full/nmeth.1237.html]. Unfortunately, it is very much built to capture Brownian motion and lacks good cost matrices for more linear movement such as axonal transport. The two files here supply this functionality. Using them I was able to achieve fully automatic quantification of mitochondrial transport in axons with R2 > 0.8 correlation compared to manual tracking.
