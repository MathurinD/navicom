#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *
VERBOSE_NAVICOM = False

nc = NaviCom("data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt")
nc.exportData("gistic")
nc.exportAnnotations()

nc.display([(('gistic', 'raw'), "size")], 'TCGA.04.1331.01')
nc.display([(('gistic', 'raw'), "glyph2_shape")], 'TCGA.04.1331.01')
# Error test
#nc.display([(('gistic', 'raw'), 'gg')])
nc.display([(('gistic', 'raw'), 'glyph1_size'), (('gistic', 'raw'), 'size')], 'TCGA.04.1331.01')

nc.display([(('gistic', 'raw'), "map_staining")])
nc.display([(('gistic', 'raw'), "map_staining"), (('gistic', 'raw'), "map_staining")]) # Second version should raise an error

nc.defineModules("data/cellcycle_v1.0.gmt")
nc.averageModule("gistic") # TODO Warning might be to fix
nc.exportData("gistic", "moduleAverage")

