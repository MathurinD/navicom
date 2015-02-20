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
        # Data, classified per analyse type
        self.processings = ["raw", "moduleAverage", "pcaComp", "geoSmooth"]
        self.data = dict() # Raw data
        self.exported_data = dict()
        for processing in self.processings:
            self.data[processing] = dict()
            self.exported_data[processing] = dict()
        self.annotations = dict() # Annotations of the samples
        self.modules = dict() # Composition of each module
        self.associated_modules = dict() # Number of modules each gene belong to
        if (fname != ""):
            self.loadData(fname)
            self.defineModules(modules_dict)
        # NaviCell export control
        self.exported_annotations = False
        self.browser_opened = False
        self.biotypes = dict()
        self.methodBiotype = {"gistic":"Discrete Copy number data", "log2CNA":"Continuous copy number data", "rna_seq_mrna":"mRNA expression data"}

    def __repr__(self):
        rpr = "NaviCom object with " + str(len(self.data)) + " types of data:\n"
        for method in self.data:
            rpr += method + ": " + self.data[method] + "\n"
        rpr += "and " + str(len(self.moduleAverage)) + " modules average:\n"
        for method in self.moduleAverage:
            rpr += method + ": " + self.moduleAverage[method] + "\n"
        return(repr)

    def loadData(self, fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt"):
        with open(fname) as file_conn:
            ff = file_conn.readlines()
            ll = 0
        # TODO Being able to import NaviCell valid files (one data type or only annotations)
        # Look for NAME or GENE in the first line and do not skip it

        while (ll < len(ff)):
            if (re.search("^M ", ff[ll])):
                # Import data
                method = re.sub("^M ", "", ff[ll].strip())
                print("Importing " + method + " data")
                samples = getLine(ff[ll+1])
                if (samples[0] == "GENE"):
                    samples = samples[1:]
                profile_data = dict()
                profile_data["samples"] = oDict()
                profile_data["genes"] = dict()
                profile_data["data"] = list()
                ii = 0
                for el in samples:
                    profile_data["samples"][el] = ii
                    ii += 1

                ll += 2
                gid = 0
                while(ll < len(ff) and not re.search("^M ", ff[ll]) and not re.search("^ANNOTATIONS", ff[ll])):
                    dl = getLine(ff[ll])
                    profile_data["data"].append(dl[1:]) # Data in a list at index gene_name
                    profile_data["genes"][dl[0]] = gid
                    gid += 1

                    ll += 1
                self.data["raw"][method] = NaviData(profile_data["data"], profile_data["genes"], profile_data["samples"])
                self.exported_data["raw"][method] = False
            elif (re.search("^ANNOTATIONS", ff[ll])):
                # Import annotations
                print("Importing Annotations")
                annotations_names = getLine(ff[ll+1])
                if (annotations_names[0] == "NAME"):
                    annotations_names = annotations_names[1:]
                annot = dict()
                annot["names"] = list()
                annot["samples"] = list()
                annot["annot"] = list()
                not_all = not "all" in annotations_names
                if (not_all):
                    annot["names"].append("all")
                    annot["annot"].append([1 for ii in range(nb_samples)])
                for name in annotations_names:
                    annot["names"].append(name)
                ll += 2
                while(ll < len(ff) and not re.search("^M ", ff[ll]) and not re.search("^ANNOTATIONS", ff[ll])):
                    al = getLine(ff[ll])
                    # Gather each annotation for this sample and add all as the first column if necessary
                    if (not_all):
                        annot["annot"].append([1] + al[1:])
                    else:
                        annot["annot"].append(al[1:])
                    annot["samples"].append(al[0])

                    ll = ll+1
                self.annotations = NaviData(annot["annot"], annot["samples"], annot["names"], dType="annotations")

    def averageModule(self, dataType):
        """
        Perform module averaging for every modules for one data type
        """
        assert(dataType in self.data["raw"].keys())
        data = self.data["raw"][dataType]
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
                    print(gene + " removed from the module")
            module_expression[module] /= non_nan
        # Put the averaging in a NaviData structure
        self.data["moduleAverage"][dataType] = NaviData(list(module_expression.values()), list(self.modules.keys()), samples)

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
        self.checkBrowser()

        for data_type in self.data:
            print(data_type) # TODO

    def exportData(self, method, processing="raw"):
        """
        Export data to NaviCell, can be processed data

        Args:
            method (str) : name of the method to export
            processing (str) : "" to export raw data, processing method to export processed data. See 'averageModule' and 'pcaComponent'
        """
        self.checkBrowser()

        if (processing in self.processings):
            if (method in self.data[processing] and method in self.methodBiotype):
                if (not self.exported_data[processing][method]):
                    self.nv.importDatatables(self.data[processing][method].makeData(self.nv.getHugoList()), method, self.methodBiotype[method])
            else:
                raise KeyError("Method " + method + " with processing " + processing + "does not exist")
        else:
            raise KeyError("Processing " + processing + " does not exist")

    def checkBrowser(self):
        if (not self.browser_opened):
            print("Launching browser...")
            self.nv.launchBrowser()
            self.browser_opened = True

    def exportAnnotations(self):
        self.checkBrowser()
        if (not self.exported_annotations):
            self.nv.sampleAnnotationImport(self.annotations.makeData())
            self.exported_annotations = True



class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name
    """

    def __init__(self, data, genes_list, samples_list, dType="data"):
        assert(len(data) == len(genes_list))
        for line in data:
            assert(len(line) == len(samples_list))
        # Initialise data and indexes
        self.data = np.array(data)
        self.rows = "genes"
        self.genes = listToDictKeys(genes_list)
        self.genes_names = list(self.genes.keys())
        self.cols = "samples"
        self.samples = listToDictKeys(samples_list)
        self.samples_names = list(self.samples.keys())
        if (dType == "annotations"):
            self.cols = "annotations"
            self.annotations = self.samples
            self.annotations_names = self.samples_names
            self.rows = "samples"
        self.dType = dType
        # For the iterator
        self.itermode = ""
        self.index = 0

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

    def __iter__(self, by="genes"):
        if (not by in ["genes", "samples"]):
            raise ValueError("'by' must be in ['genes', 'samples']")
        self.iter_mode = by
        self.index = 0
        return(self)

    def __next__(self):
        try:
            if (self.iter_mode == "genes"):
                key = list(self.genes.keys())[self.index]
                result = [key] + list(self.data[self.index,:])
            elif (self.iter_mode == "samples"):
                key = list(self.samples.keys())[self.index]
                result = [key] + list(self.data[:,self.index])
        except IndexError:
            raise StopIteration
        self.index += 1
        return(result)

    def __repr__(self):
        rpr = "NaviData array with " + str(len(self.genes)) + " " + self.rows + " and "
        rpr += str(len(self.samples)) + " " + self.cols
        return(rpr)

    def makeData(self, hugo_map=""):
        """ Builds a string suitable for NaviCell Web Service from a python matrix of gene/sample values or a NaviCom object.

        Matrix format:
        - first line is: GENE word followed by a tab separated list of sample names,
        - each line begins with an gene name and must be followed by a tab separated list of gene/sample values.

        Eliminates genes not present in the NaviCell map.
        """

        # Change header whether we have annotations or real self
        if (self.dType == "annotations"):
            ret = "NAME\t" + '\t'.join(self.annotations_names) + "\n"
            for line in self:
                ret += buildLine(line)
        elif (self.dType == "data"):
            ret = "GENE\t" + '\t'.join(self.samples_names) + "\n"
            for line in self:
                if (line[0] in hugo_map or hugo_map == ""):
                    ret += buildLine(line)
        else:
            raise ValueError("Data type must be 'data' or 'annotations'")

        return("@DATA\n" + ret)
            

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

    def __iter__(self):
        return(self.data.__iter__())

    def __repr__(self):
        rpr = "NaviData slice with " + str(len(self.ids)) + " elements"
        return(rpr)

    def __add__(self, nvslice):
        assert(isinstance(nvslice, NaviSlice))
        #assert(self.ids == nvslice.ids)
        return(NaviSlice(self.data + nvslice.data, self.ids))

def listToDictKeys(ilist):
    if (isinstance(ilist, dict)):
        return(oDict(ilist))
    res = oDict()
    ii = 0
    for key in ilist:
        res[str(key)] = ii
        ii += 1
    return(res)

def buildLine(line, sep="\t"):
    ret = ""
    for el in line:
        ret += str(el) + sep
    return(ret[:-1] + "\n")




