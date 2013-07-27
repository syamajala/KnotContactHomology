#!/usr/bin/sage -python
from braids_parallel import *

print "Make sure you run sage -sh; ipcluster start --profile=sage;"

setup_view()
braid = Braid(satellite([1, 1, 1]))
print "Braid: " + str(satellite([1, 1, 1]))
braid.homology()
