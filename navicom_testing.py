#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *
VERBOSE_NAVICOM = False

nc = NaviCom()
nc.defineModules("data/cellcycle_v1.0.gmt")
nc.averageModule("gistic") # TODO Warning might be to fix

nc.exportData("gistic")
nc.exportData("gistic", "moduleAverage")
nc.exported_annotations = False
nc.exportAnnotations()

nc.display([(('gistic', 'raw'), "map_staining")])

