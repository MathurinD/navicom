#########################################################################################
# navicomgui.py
# By Mathurin Dorel
# A GUI to download cBioPortal data and use them in NaviCom
#########################################################################################

import os
import tkinter as tk
import tkinter.filedialog as fd
import tkinter.colorchooser as cc

from navicom import *
from selectionlist import *

main = tk.Tk()
main.geometry("{}x{}+100+100".format(main.winfo_screenwidth(), main.winfo_screenheight()))
main.title("NaviCom")

options = ["ption 1", "null", "value", "because"]
sample_selection = SelectionList(main, options)
map_selection = SelectionList(main, ["sample1", "sample2"])

#main.iconify()
main.mainloop()

