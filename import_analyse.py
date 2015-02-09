# Import data from a file and create analyzed sets to export
from curie.navicell import *
import numpy as np
import re
from collections import OrderedDict as oDict

def getLine(ll, split_char="\t"):
    ll = re.sub('\"', '', ll.strip())
    ll = re.sub('\tNA', '\tNaN', ll) # Python float can convert "NaN" -> nan but not "NA"
    ll = ll.split(split_char)
    for ii in ll:
        try:
            ll[ii] = float(ll[ii])
        except:
            pass
    return(ll)

class NaviCom():

    def __init__(self, fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt"):
        self.data = dict()
        self.annotations = dict()
        if (fname != ""):
            try:
                self.loadData(fname)
            except:
                print("Data could not be loaded, empty object returned")

    def loadData(self, fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt"):
        with open(fname) as file_conn:
            ff = file_conn.readlines()
            ll = 0

        all_data = dict()
        annot_dict = dict()
        while (ll < len(ff)):
            if (re.search("^M ", ff[ll])):
                # Import data
                method = re.sub("^M ", "", ff[ll].strip())
                print("Importing " + method)
                samples = getLine(ff[ll+1])
                profile_data = dict()
                profile_data["samples"] = oDict()
                ii = 0
                for el in samples:
                    profile_data["samples"][el] = ii
                    ii += 1

                ll += 2
                while(ll < len(ff) and not re.search("^M ", ff[ll]) and not re.search("^ANNOTATIONS", ff[ll])):
                    dl = getLine(ff[ll])
                    profile_data[dl[1]] = dl[:-1] # Data in a list at index gene_name

                    ll += 1
                all_data[method] = profile_data
            elif (re.search("^ANNOTATIONS", ff[ll])):
                # Import annotations
                print("Importing Annotations")
                annotations_names = getLine(ff[ll+1])
                annot_dict["annotations"] = oDict()
                ii = 0
                for annot in annotations_names:
                    annot_dict["annotations"][annot] = ii
                    annot_dict[annot] = list()
                    ii += 1

                ll += 2
                while(ll < len(ff) and not re.search("^M ", ff[ll]) and not re.search("^ANNOTATIONS", ff[ll])):
                    al = getLine(ff[ll])
                    # Gather each annotation for this sample
                    spl = al[1]
                    for name in annotations_names:
                        annot = al[annot_dict["annotations"][name]]
                        annot_dict[name].append(annot)
                    ll = ll+1
                if (not "all" in annotations_names):
                    annot_dict["annotations"].append["all"]
                    annot_dict["all"].append(1)
        self.data = all_data
        self.annotations = annot_dict



