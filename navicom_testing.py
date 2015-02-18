#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *

nc = NaviCom()
nc.defineModules("data/cellcycle_v1.0.gmt")
nc.averageModule("gistic")

print(nc.moduleAverage)
nc.data['gistic']['PURA', 'AKT1']
#nc.data['gistic']['PURA', 'truc'] # Error generating test

