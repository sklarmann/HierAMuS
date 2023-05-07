# Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

import os, sys
import gmsh


gmsh.initialize()


pathname = os.path.dirname(sys.argv[0])

fname = os.path.join(pathname,'testgeo.geo')

gmsh.open(fname)



gmsh.fltk.run()
