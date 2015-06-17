#########################################################################################
# navicomgui.py
# By Mathurin Dorel
# A GUI to download cBioPortal data and use them in NaviCom
#########################################################################################

import os
import re
import tkinter as tk
import tkinter.filedialog as fd
import tkinter.colorchooser as cc

from navicom import *
from selectionlist import *

def skipRHeader(rexec):
    """ Return the R result stream after R headers and the command line """
    ll = rexec.readline().strip()
    while (not re.search("library", ll)):
        ll = rexec.readline().strip()

def listStudies():
    """ Retrieve the list of cBioPortal studies that are available """
    rprogram = """library(cBioFetchR);options("max.print"=10000)
    conn = cBioConnect()
    listStudies(conn)
    """
    rprogram = re.sub("\n", ";", rprogram.strip())
    studies = dict()
    studies["id"] = list()
    studies["name"] = list()
    # Use R to get the study
    with os.popen("R -e '" + rprogram + "'") as rexec:
        skipRHeader(rexec)
        ll = rexec.readline().strip()
        ll = rexec.readline().strip()
        print(ll)
        ll = re.sub(" +([^ ]+) +(.+)$", "\t\\1\t\\2", rexec.readline().strip())
        while (re.match("^\d", ll)):
            sl = ll.split("\t")
            print(sl)
            studies["id"].append(sl[1])
            studies["name"].append(sl[2])
            print(ll)
            ll = re.sub(" +([^ ]+) +(.+)$", "\t\\1\t\\2", rexec.readline().strip())
    return(studies)

main = tk.Tk()
main.geometry("{}x{}+100+100".format(main.winfo_screenwidth(), main.winfo_screenheight()))
main.title("NaviCom")

studies = listStudies()
study_selection = SelectionList(main, "Study :", studies["name"], studies["id"])

map_selection = SelectionList(main, "Map selection", ["sample1", "sample2"])

method_selection = SelectionList(main, "Selection of the display method", ["Complete", "Methylation", "Omics"])

#main.iconify()
main.mainloop()

