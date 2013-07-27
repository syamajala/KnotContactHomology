#!/usr/bin/sage -python
from braids import *

braid = Braid([1])
print "Braid: [1]"
braid.homology()

braid = Braid(satellite([1]))
print "Braid: " + str(satellite([1]))
braid.homology()

braid = Braid([1, 1, 1])
print "Braid: [1, 1, 1]"
braid.homology()

