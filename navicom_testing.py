#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *
import time
VERBOSE_NAVICOM = False

nc = NaviCom("data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt")
nc.exportAnnotations() # TODO move annotations export check to data export function
nc.dataAvailable()

nc.displayMethylome(['TCGA.04.1331.01'], "raw", "mRNA", "size")

nc.display([('log2CNA', 'barplot')], ['OS_STATUS: NA; SEQUENCED: NA'])

nc.resetDisplay()

nc.display([('log2CNA', 'barplot')], 'OS_STATUS: NA')
nc.display([('log2CNA', 'barplot')], ['OS_STATUS; SEQUENCED', 'all_groups'])

nc.display([('log2CNA', 'heatmap')], 'OS_STATUS: NA')
nc.display([('log2CNA', 'heatmap'), ('gistic', 'heatmap')], ['OS_STATUS', 'all_groups'])

nc.exportData("gistic")

nc.display([('log2CNA', 'barplot')], 'TCGA.04.1331.01')

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

