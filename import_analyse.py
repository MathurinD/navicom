# Import data from a file and create analyzed sets to export
from curie.navicell import *
import numpy as np
import re
from collections import OrderedDict as oDict
from pprint import pprint

DEBUG_NAVICOM = True
VERBOSE_NAVICOM = True

def getLine(ll, split_char="\t"):
    ll = re.sub('\"', '', ll.strip())
    ll = re.sub('\tNA', '\tNaN', ll) # Python float can convert "NaN" -> nan but not "NA" -> nan
    ll = ll.split(split_char)
    for ii in range(len(ll)):
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
        self.data = dict() # Raw data
        self.annotations = dict() # Annotations of the samples
        self.modules = dict() # Composition of each module
        self.associated_modules = dict() # Number of modules each gene belong to
        if (fname != ""):
            self.loadData(fname)
            self.defineModules(modules_dict)
        # Analysed structures, indexed by data type
        self.moduleAverage = dict()
        #self.pcaAnalyse = dict()

    def __repr__(self):
        rpr = "NaviCom object with " + len(self.data) + " types of data:\n"
        for method in self.data:
            rpr += method + ": " + self.data[method] + "\n"
        return(repr)

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
                print("Importing " + method + " data")
                samples = getLine(ff[ll+1])
                profile_data = dict()
                profile_data["samples"] = dict()
                ii = 0
                for el in samples:
                    profile_data["samples"][el] = ii
                    ii += 1

                ll += 2
                profile_data["genes"] = dict()
                profile_data["data"] = list()
                gid = 0
                while(ll < len(ff) and not re.search("^M ", ff[ll]) and not re.search("^ANNOTATIONS", ff[ll])):
                    dl = getLine(ff[ll])
                    profile_data["data"].append(dl[1:]) # Data in a list at index gene_name
                    profile_data["genes"][dl[0]] = gid
                    gid += 1

                    ll += 1
                all_data[method] = NaviData(profile_data["data"], profile_data["genes"], profile_data["samples"])
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
                    for line in ff.readlines():
                        ll = line.strip().split("\t")
                        module_name = ll[0]
                        self.modules[module_name] = list()
                        if (ll[1] != "na"):
                            self.modules[module_name] += ll[1]
                        self.modules[module_name] += ll[2:]
                        # Count the number of modules each gene belong to
                        for gene in self.modules[module_name]:
                            try:
                                self.associated_modules[gene] += 1
                            except KeyError:
                                self.associated_modules[gene] = 1
        """
        # Only keep genes with data
        for module in self.modules:
            keep = list()
            for gene in self.modules[module]:
                if (not gene in self.genes_list):
                    self.modules[module].remove(gene)
        """

    def display(self):
        """
        Display the data on a NaviCell map
        """
        self.launchBrowser()

        for data_type in self.data:
            print(data_type) # TODO

    def averageModule(self, dataType):
        """
        Perform module averaging for every modules for one data type
        """
        assert(dataType in self.data.keys())
        data = self.data[dataType]
        samples = list(data.samples.keys())
        module_expression = dict()
        for module in self.modules:
            module_expression[module] = [0 for sample in samples]
            non_nan = np.array([0 for sample in samples])
            no_data = list()
            for gene in self.modules[module]:
                try:
                    not_nan = [int(not np.isnan(dd)) for dd in data[gene].data]
                    non_nan += not_nan
                    module_expression[module] += data[gene].data * not_nan
                except IndexError:
                    no_data += [gene]
            for gene in no_data:
                self.modules[module].remove(gene)
                if (VERBOSE_NAVICOM):
                    print(gene + " removed")
            module_expression[module] /= non_nan
        # Put the averaging in a NaviData structure
        self.moduleAverage[dataType] = NaviData(list(module_expression.values()), list(self.modules.keys()), samples)



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
                raise IndexError("'" + index + "' is neither a gene or sample name")
        elif (isinstance(index, list) or isinstance(index, tuple)):
            result = list()
            for ii in index:
                result.append(self[ii].data)
            result = np.array(result)
            if (index[0] in self.genes):
                assert np.all([idx in self.genes for idx in index]), "Not all index are gene names"
                genes = self.genes.copy()
                for gene in self.genes:
                    if (not gene in index):
                        genes.pop(gene)
                return( NaviData(result, genes, self.samples) )
            elif (index[0] in self.samples):
                assert np.all([idx in self.samples for idx in index]), "Not all index are sample names"
                samples = self.samples.copy()
                for sample in self.samples:
                    if(not sample in index):
                        samples.pop(sample)
                return( NaviData(result.transpose(), self.genes, samples) )

    def __repr__(self):
        rpr = "NaviData array with " + str(len(self.genes)) + " genes and "
        rpr += str(len(self.samples)) + " samples"
        return(rpr)
            

class NaviSlice():
    """
    A slice from a NaviData array
    """
    def __init__(self, data, index_dict):
        if (isinstance(data, np.ndarray)):
            self.data = data
        elif (isinstance(data, list)):
            self.data = np.array(data)
        else:
            raise TypeError
        if (isinstance(index_dict, dict)):
            self.ids = index_dict
        elif (isinstance(index_dict, list)):
            self.ids = listToDictKeys(index_dict)
        else:
            raise TypeError

    def __getitem__(self, index):
        if (isinstance(index, int)):
            return(self.data[index])
        elif (isinstance(index, str)):
            return(self.data[self.ids[index]])
        elif (isinstance(index, list)):
            for ii in index:
                return(self[ii])

    def __repr__(self):
        rpr = "NaviData slice with " + str(len(self.ids)) + " elements"
        return(rpr)

    def __add__(self, nvslice):
        assert(isinstance(nvslice, NaviSlice))
        #assert(self.ids == nvslice.ids)
        return(NaviSlice(self.data + nvslice.data, self.ids))

def listToDictKeys(ilist):
    if (isinstance(ilist, dict)):
        return(ilist)
    res = dict()
    ii = 0
    for key in ilist:
        res[str(key)] = ii
        ii += 1
    return(res)




