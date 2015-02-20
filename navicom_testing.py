#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *
nv.VERBOSE_NAVICOM = False

nc = NaviCom()
nc.defineModules("data/cellcycle_v1.0.gmt")
nc.averageModule("gistic")

nc.checkBrowser()
nc.nv.importDatatables(nc.nv.makeDataFromFile("data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011_gistic.tsv"), "gistic", "Discrete Copy number data")
nc.exported_annotations = False
nc.exportAnnotations()

