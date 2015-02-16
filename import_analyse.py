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

    def __init__(self, fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt", map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php', modules_dict=""):
        # Build options for the navicell connexion
        options = Options()
        options.map_url = map_url
        idx = options.map_url.find('/navicell/')
        options.proxy_url = options.map_url[0:idx] + '/cgi-bin/nv_proxy.php'
        options.browser_command = "firefox %s" # TODO Add user control
        self.nv = NaviCell(options)
        # Data
        self.data = dict()
        self.annotations = dict()
        if (fname != ""):
            try:
                self.loadData(fname)
            except:
                print("Data could not be loaded, empty object returned")
        self.defineModules(modules_dict)

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
                samples = getLine(ff[ll+1])[1:]
                profile_data = dict()
                profile_data["samples"] = oDict()
                ii = 0
                for el in samples:
                    profile_data["samples"][el] = ii
                    ii += 1

                ll += 2
                profile_data["genes"] = dict()
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

    def defineModules(self, modules_dict=""):
        """
        Defines the modules to use and which module each gene belongs to

        Args:
            modules_dict : Either a dict indexed by module name or a file name with the description of each module (.gmt file)
        """
        self.modules = dict()
        if (isinstance(modules_dict, dict)):
            self.modules = modules_dict # TODO add control that the genes are included
        elif (isinstance(modules_dict, str)):
            if (modules_dict != ""):
                # TODO get the list of modules from the file
                with open(modules_dict) as ff:
                    for line in ff.readLines():
                        ll = line.strip().split("\t")
                        module_name = ll[0]
                        self.modules[module_name] = list()
                        if (ll[1] != "na"):
                            self.modules[module_name].append(ll[1])
                        self.modules[module_name].append(ll[2:])
        # Only keep genes with data
        for module in self.modules:
            keep = list()
            for gene in self.modules[module]:
                if (gene in self.genes_list):
                    self.modules[module].remove(gene)

    def display(self):
        """
        Display the data on a NaviCell map
        """
        self.launchBrowser()

        for data_type in self.data:
            print(data_type) # TODO

    #def 

class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name
    """

    def __init__(self, data, genes_list, samples_list):
        assert(len(data) == len(genes_list))
        for line in data:
            assert(len(line) == len(samples_list))
        # Initialise data and indexes
        self.data = np.array(data)
        self.genes = listToDictKeys(genes_list)
        self.samples = listToDictKeys(samples_list)

    def __getitem__(self, index):
        if (isinstance(index, int)):
            return(self.data[index])
        elif (isinstance(index, str)):
            if (index in self.genes):
                return( NaviSlice(self.data[self.genes[index],:], self.samples) )
            elif (index in self.samples):
                return( NaviSlice(self.data[:,self.samples[index]], self.genes) )
            else:
                raise IndexError

    def __repr__(self):
        rpr = "NaviData array with " + str(len(self.genes)) + " genes, and "
        rpr += str(len(self.samples)) + " samples"
        return(rpr)
            

class NaviSlice():
    """
    A slice from a NaviData array
    """
    def __init__(self, data, index_dict):
        self.data = data
        self.ids = index_dict

    def __getitem__(self, index):
        if (isinstance(index, int)):
            return(self.data[index])
        elif (isinstance(index, str)):
            return(self.data[self.ids[index]])

    def __repr__(self):
        rpr = "NaviData slice with " + str(len(self.ids)) + " elements"
        return(rpr)

def listToDictKeys(ilist):
    if (isinstance(ilist, dict)):
        return(ilist)
    res = dict()
    ii = 0
    for key in ilist:
        res[str(key)] = ii
        ii += 1
    return(res)




