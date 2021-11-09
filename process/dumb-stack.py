from glob import glob
import numpy as np
from astropy.nddata import CCDData
import ccdproc

files = glob('stacks/c2014un271/*med.fits')
images = [CCDData.read(f, unit='adu') for f in files]
combiner = ccdproc.Combiner(images)
combiner.sigma_clipping()
avg = combiner.average_combine()
med = combiner.median_combine()
