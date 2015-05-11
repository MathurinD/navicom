#!/usr/bin/python3
#-*-coding:utf-8-*-

from import_analyse import *

nc = NaviCom(map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php', fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt")
nc.listData()
nc.listAnnotations()

nc.colorsOverlay("mrna_median", "log2CNA", processing="raw")
nc.listData()
nc.saveData( "mrna_median_log2CNA", "colors")

nc.saveAllData()
nc.loadData("data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011_gistic.tsv")
nc.loadData("Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.ncc")

nc.displayTranscriptome('log2CNA', 'OS_STATUS: LIVING', "barplot", 'quantiles')

dd=nc.generateDistributionData(nc.getDataName('log2CNA'), 'OS_STATUS: LIVING')
nc.distData[dd[0]].exportToNaviCell(nc.nv, TYPES_BIOTYPE['mRNA'], dd[0])

nc.exportAnnotations()

nc.displayTranscriptome('log2CNA', 'OS_STATUS: NA', "barplot", 'TCGA.04.1331.01')

nc.defineModules("data/cellcycle_v1.0.gmt")
nc.averageModule("gistic") # TODO Warning might be to fix

nc.display([('log2CNA', 'barplot'), ('gistic', 'shape2', 'TCGA.04.1331.01'), ('log2CNA', 'size2', 'TCGA.04.1331.01')], ['OS_STATUS: NA; SEQUENCED: NA'])

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

