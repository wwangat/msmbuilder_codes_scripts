#calculate the pearson correlation (the one deduct mean) of tica w.r.t. the known reaction coordinates
from pydoc import help
from scipy.stats.stats import pearsonr
import numpy

for line in file('pearson.file'):
    data = numpy.loadtxt(line.strip())
    print line.strip(), pearsonr(data[:, 0], data[:, 1])
