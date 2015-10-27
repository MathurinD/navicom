################################################################################
# navidata.py
# By Mathurin Dorel
# Classes to store data an annotations in a convenient format for NaviCell export
################################################################################

import numpy as np
import re
from time import sleep
from collections import OrderedDict as oDict
from pprint import pprint
import math
from warnings import warn
from .displayConfig import *

DEBUG_NAVIDATA = True

# Constants related to NaviCell
MAX_GLYPHS = 5
GLYPH_TYPES = ["color", "size", "shape"]
# Identify the categories of cBioportal data as (aliases, biotype)
# Also add some generic designations for each type of biological data
# The first in the list are prefered for display when several are available
TYPES_SPEC = dict()
TYPES_SPEC["mRNA"] = (["mrna", "rna_seq_v2_mrna", "rna_seq_mrna", "rna_seq_rna", "zscores", "mrna_median", "mrna_median_zscores", "rna_seq_mrna_median_zscores", "rna_seq_v2_mrna_median_zscores", "mrna_U133", "mrna_U133_zscores", "mrna_zbynorm", "mrna_outliers", "mrna_znormal", "mrna_outlier", "mrna_merged_median_zscores"], "mRNA expression data", "none")
TYPES_SPEC["dCNA"] = (["gistic", "cna", "cna_rae", "cna_consensus", "snp-fasst2"], "Discrete Copy number data", "none")
TYPES_SPEC["cCNA"] = (["log2cna"], "Continuous copy number data", "none")
TYPES_SPEC["methylation"] = (["methylation", "methylation_hm27", "methylation_hm450"], "mRNA expression data", ("ff7000", "diamond"))
TYPES_SPEC["protein"] = (["protein_level", "rppa_protein_level", "proteomics", "rppa", "rppa_zscores"], "Protein Expression Data", ("ffff00", "circle"))
TYPES_SPEC["miRNA"] = (["mirna", "mirna_median_zscores"], "microRNA expression data", ("8800ff", "hexagon"))
TYPES_SPEC["mutations"] = (["mutations"], "Continuous copy number data", ("0000ff", "triangle"))
#TYPES_SPEC["mutations"] = (["mutations"], "Mutations")
TYPES_SPEC["unknown"] = (["unknown"], "mRNA expression data", "none") # If the type of data cannot be identified, consider continuous data by default
METHODS_TYPE = dict()
TYPES_BIOTYPE = dict()
TYPE_CONFIG = dict()
for bt in TYPES_SPEC:
    TYPES_BIOTYPE[bt] = TYPES_SPEC[bt][1]
    TYPE_CONFIG[bt] = TYPES_SPEC[bt][2]
    for method in TYPES_SPEC[bt][0]:
        METHODS_TYPE[method.lower()] = bt

# Inform what the processing does to the biotype, if -> it changes only some biotypes, if "something" it turns everything to something, see exportData and saveData
PROCESSINGS = ["raw", "moduleAverage", "pcaComp", "geoSmooth", "distribution", "colors", "textMutations"]
PROCESSINGS_BIOTYPE = {"moduleAverage":"Discrete->Continuous", "pcaComp":"Color", "geoSmooth":"Discrete->Continuous", "textMutations":"Mutation data", "distribution":"Continuous Copy number data"}
DISCRETE_BIOTYPES = ["Mutation data", "Discrete Copy number data"]
CONTINOUS_BIOTYPES = ["mRNA expression data", "microRNA expression data", "Protein Expression Data", "Continuous copy number data"]

class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name.
    Args :
        data (list or array) : Values of the data to insert in the NaviData object. Must be convertible into a numpy array.
        rows_list (list) : names of the rows (samples names)
        columns_list (list) : names of the columns (genes names)
        processing (str) : name of the computer processing applied to the data
        method (str) : name of the experimental method used to get the original ("raw") data
        dType (str) : "data" or "annotations", whether the NaviData object contains datas or annotations (Note : this should be left to default, this is used by NaviAnnotations to change some internal variables)
        display_config (GlyphConfig): An optionnal configuration if those data are to be used as glyphs, to set a constant color and glyph symbol
    """
    def __init__(self, data, rows_list, columns_list, method="unknown", processing="raw", dType="data", display_config=""):
        assert(len(data) == len(rows_list))
        for line in data:
            assert len(line) == len(columns_list), "Incorrect length of line : " + str(line) + ", length = " + str(len(line)) + ", expected " + str(len(columns_list))
        # Informations on the data TODO
        self.processing = processing
        self.method = method
        self.biotype = getBiotype(method, processing)
        # Initialise data and indexes
        ## Raw datatable as a numpy array
        self.data = np.array(data)
        ## Names of the rows
        self._rows = listToDictKeys(rows_list)
        self.rows_names = list(self._rows.keys())
        ## Names of the columns
        self._columns = listToDictKeys(columns_list)
        self.columns_names = list(self._columns.keys())
        if (dType == "annotations"):
            self.inColumns = "annotations"
            self._annotations = self._columns
            self.annotations_names = self.columns_names
            self.inRows = "samples"
            self._samples = self._rows
            self.samples_names = self.rows_names
        elif (dType == "data"):
            self.inRows = "genes"
            self._genes = self._rows
            self.genes_names = self.rows_names
            self.inColumns = "samples"
            self._samples = self._columns
            self.samples_names = self.columns_names
        elif (dType == "old_annotations"):
            self.inRows = "annotations"
            self.inColumns = "samples"
        else:
            raise ValueError("dType '" + dType + "' is not valid")
        self.dType = dType
        # For the iterator
        self.itermode = ""
        self.index = 0
        if (isinstance(display_config, GlyphConfig)):
            self.display_config = display_config
        elif (display_config == "gradient"):
            self.display_config = display_config
        elif (display_config == ""):
            config = TYPE_CONFIG[getDataType(method, processing)]
            if (config == "none"):
                self.display_config = "gradient"
            else:
                self.display_config = GlyphConfig(config[0], config[1]) # Use predefined configs
        else:
            raise ValueError("Invalid display configuration: " + display_config)

    def __getitem__(self, index):
        if (isinstance(index, int)):
            return(self.data[index])
        elif (isinstance(index, str)):
            if (index in self._rows):
                return( NaviSlice(self.data[self._rows[index],:], self._columns) )
            elif (index in self._columns):
                return( NaviSlice(self.data[:,self._columns[index]], self._rows) )
            else:
                raise IndexError("'" + index + "' is neither a gene or sample name")
        elif (isinstance(index, list) or isinstance(index, tuple)):
            result = list()
            for ii in index:
                result.append(self[ii].data)
            result = np.array(result)
            if (index[0] in self._rows):
                assert np.all([idx in self._rows for idx in index]), "Not all index are gene names"
                genes = self._rows.copy()
                for gene in self._rows:
                    if (not gene in index):
                        genes.pop(gene)
                return( NaviData(result, genes, self._columns) )
            elif (index[0] in self._columns):
                assert np.all([idx in self._columns for idx in index]), "Not all index are sample names"
                samples = self._columns.copy()
                for sample in self._columns:
                    if(not sample in index):
                        samples.pop(sample)
                return( NaviData(result.transpose(), self._rows, samples) )

    def __iter__(self, by="genes"):
        if (not by in ["genes", "samples"]):
            raise ValueError("'by' must be in ['genes', 'samples']")
        self.iter_mode = by
        self.index = 0
        return(self)

    def __next__(self):
        try:
            if (self.iter_mode == "genes"):
                key = list(self._rows.keys())[self.index]
                result = [key] + list(self.data[self.index,:])
            elif (self.iter_mode == "samples"):
                key = list(self._columns.keys())[self.index]
                result = [key] + list(self.data[:,self.index])
        except IndexError:
            raise StopIteration
        self.index += 1
        return(result)

    def __repr__(self):
        rpr = "NaviData array with " + str(len(self._rows)) + " " + self.inRows + " and "
        rpr += str(len(self._columns)) + " " + self.inColumns
        return(rpr)

    def _makeData(self, hugo_map=""):
        """ Builds a string suitable for NaviCell Web Service from a python matrix of gene/sample values or a NaviCom object.

        Matrix format:
        - first line is: GENE word followed by a tab separated list of sample names,
        - each line begins with an gene name and must be followed by a tab separated list of gene/sample values.

        Remove genes not present in hugo_map if provided.
        """

        # Change header whether we have annotations or real self
        if (self.dType == "annotations"):
            ret = "NAME\t" + '\t'.join(self.annotations_names) + "\n"
            for line in self:
                ret += buildLine(line)
        elif (self.dType == "data"):
            ret = "GENE\t" + '\t'.join(self.columns_names) + "\n"
            for line in self:
                if (line[0] in hugo_map or hugo_map == ""):
                    ret += buildLine(line)
        else:
            raise ValueError("Data type must be 'data' or 'annotations'")

        return("@DATA\n" + ret)

    def exportToNaviCell(self, nv, biotype, dataName):
        """
        Export the datatable to a NaviCell map
        
        Args:
            nv (NaviCell): a NaviCell communication object
            biotype (str): biotype of the datatable in NaviCell
            dataName (str): name of the datatable in NaviCell
        """
        nv.importDatatables(self._makeData(nv.getHugoList()), dataName, biotype) # Remove name once in self
        #nv.importDatatables(self._makeData(nv.getHugoList()), dataName, self.biotype)

    def saveData(self, baseName="", mode="w", fullName=False):
        """
        Save the NaviData datas in a file that can be used in NaviCell, but can also be loaded as NaviData

        Args:
            baseName (str): first part of the name of the file where the data will be written. The processing and method are added to file name if fullName is False.
            mode (str): How the file should be opened ('a' or 'w')
            fullName (bool): Whether the baseName should be used alone (True) or with the extra information of the method and processing.
        """
        assert mode in ["a", "w"], ValueError("Cannot open the file with this mode to save data")
        # Build filename
        if (self.dType == "data"):
            fname = baseName + "[" + self.processing + "]" + "[" + self.method + "].ncd"
        elif (self.dType == "annotations"):
            fname = baseName + "_annotations.nca"
        # Save data
        if (mode == "a" or fullName):
            fname = baseName
        with open(fname, mode) as ff:
            # First line, and first word in the 
            if (self.dType == "data"):
                if (mode == "a"):
                    ff.write("M "+self.method+"\t"+self.processing+"\n")
                ff.write("GENE")
            elif (self.dType == "annotations" or self.dType == "old_annots"):
                if (mode == "a"):
                    ff.write("ANNOTATIONS "+self.method+"\t"+self.processing+"\n")
                ff.write("NAME")
            for nn in self.columns_names:
                ff.write("\t" + nn)
            ff.write("\n")
            for ii in range(len(self.rows_names)):
                ff.write(self.rows_names[ii])
                for sv in self.data[ii]:
                    ff.write("\t" + str(sv))
                ff.write("\n")
            
MAX_GROUPS = 7
class NaviAnnotations(NaviData):
    """
    Enhance NaviData to contain annotations and associate annotations values with samples. Also reduce continuous data with to many levels to a limited number of interval levels.
    """

    def __init__(self, data, rows_list, columns_list, dType="annotations"):
        NaviData.__init__(self, data, rows_list, columns_list, dType=dType)
        # Get the values for each annotations and which samples are associated to a value
        self.categoriesPerAnnotation = dict()
        self._samplesPerCategory = dict()
        modified_data = list()
        modified_annot = list()
        for annot in self._annotations:
            self.categoriesPerAnnotation[annot] = list()
            self._samplesPerCategory[annot] = dict()
            reduced = False
            # Reduce the annotations set if they are continuous integers with too many values
            if (len(np.unique(self[annot].data)) > 1.5 * MAX_GROUPS):
                try:
                    values = [float(value) for value in np.unique(self[annot].data)] # ValueError if not floats
                    # Define the new annotations
                    min_value = min(values)
                    max_value = max(values)
                    step = (max_value-min_value)/MAX_GROUPS
                    new_values = [signif(cat,4) for cat in np.arange(min_value, max_value+step/2, 1.01*step)]
                    new_categories = [str(new_values[icat])+"<X<"+str(new_values[icat+1]) for icat in range(len(new_values)-1)] + ["NaN"]
                    self.categoriesPerAnnotation[annot] = new_categories
                    # Attribute the new annotations to samples
                    modified_data.append(self[annot].data.copy())
                    modified_annot.append(annot)
                    for cat in new_categories:
                        self._samplesPerCategory[annot][cat] = list()
                    for sample in self._samples:
                        annot_value = float(self[annot][sample])
                        if (np.isnan(annot_value)):
                            self._samplesPerCategory[annot]["NaN"].append(sample)
                        else:
                            icat = 0
                            while (signif(annot_value) >= new_values[icat+1]):
                                icat += 1
                            self._samplesPerCategory[annot][new_categories[icat]].append(sample)
                            self[annot][sample] = new_categories[icat]
                    reduced = True
                except ValueError: 
                    reduced = False
            # Non numeric values or with few enough levels are not modified and simply indexed
            if (not reduced):
                for sample in self._samples:
                    annot_value = self[annot][sample]
                    if (not annot_value in self.categoriesPerAnnotation[annot]):
                        self.categoriesPerAnnotation[annot].append(annot_value)
                        self._samplesPerCategory[annot][annot_value] = [sample]
                    else:
                        self._samplesPerCategory[annot][annot_value].append(sample)

        if (DEBUG_NAVIDATA):
            print("Discretised annotations:" + str(modified_annot))
        if (len(modified_annot) > 0):
            self.old_annots = NaviData(modified_data, modified_annot, self._rows, dType="old_annotations")
        else:
            self.old_annots = list()

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
            ret_value = list()
            ret_name = list()
            for ii in index:
                if (ii in self.ids):
                    ret_name.append(ii)
                else:
                    ret_name.append(self.ids[ii])
                ret_value.append(self[ii])
            return(NaviSlice(ret_value, ret_name))

    def __setitem__(self, index, value):
        if (isinstance(index, int)):
            self.data[index] = value
        elif (isinstance(index, str)):
            self.data[self.ids[index]] = value
        elif (isinstance(index, list)):
            if (isinstance(value, list) and len(value) == len(index)):
                for ii in range(len(index)):
                    self.data[index[ii]] = value[ii]
            else:
                raise ValueError("Length of indexes to change and new values differ")

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
    """
    Convert a list to an ordered dict where each item is a key and gives the id of the item in the list
    """
    if (isinstance(ilist, dict)):
        return(oDict(ilist))
    res = oDict()
    ii = 0
    for key in ilist:
        res[validateDimNames(str(key))] = ii
        ii += 1
    return(res)

def validateDimNames(dimname):
    """
        Make sure the names of the rows and columns are valid to be used as javascript selectors
    """
    out = re.sub("#|\.|@|\+", "_", dimname)
    return( out )

def buildLine(line, sep="\t"):
    ret = ""
    for el in line:
        ret += str(el) + sep
    return(re.sub("NaN|nan", "NA", ret[:-1]) + "\n")

def parseFilename(fname):
    """
    Extract data name, method and processing from the name of the file
    """
    fname = re.sub(".*/([^/]+)", "\\1", fname)
    if (re.search("\.ncd$|\.nca$", fname)):
        names = re.split("\[", re.sub(".ncd$|.nca$", "", fname))
        dname = re.sub("_$", "", names[0])
        processing = re.sub("\]|_$|^_", "", names[1])
        method = re.sub("\[|_$|^_", "", names[2])
    elif (re.search("\.tsv$", fname)):
        processing = "raw" # No processing is done by the R extractor
        names = re.split("_", re.sub("\.tsv", "", fname))
        ll = len(names)
        while (not "_".join(names[-ll:]) in METHODS_TYPE and ll >= 1):
            ll-=1
        if (ll >= 1):
            method = "_".join(names[-ll:])
            dname = "_".join(names[:-ll])
    else:
        return (re.sub("\.[^.]*$", "", fname), "raw", "unknown")
    return (dname, processing, method)

def pca_cov(data):
    """
    Run Principal Component Analysis using the covariance matrix and returns the eigen-vectors sorted by decreasing eigenvalues
    """
    cov_mat = np.cov(data)

    # Compute eigenvectors and eigenvalues, linalg returns eigenvectors of norm 1
    eig_val, eig_vec = np.linalg.eig(cov_mat)
    # Make a list of value-vector tuples to sort by eigenvalue
    eig_pairs = [(np.abs(eig_val[ii]), eig_vec[:,ii]) for ii in range(len(eig_val))]
    eig_pairs.sort()
    eig_pairs.reverse()

    # Create the ordered eigenvectors matrix
    eigen_matrix = np.matrix([eig_pairs[1] for ii in range(len(eig_pairs))])
    return (eigen_matrix)

def signif(x, n=3):
    """
    Keep n significant numbers
    """
    if (x==0):
        return 0.
    return(round(x, -int(math.log10(np.abs(x)))+(n-1) ))

def getBiotype(method, processing="raw"):
    """
        Get the biotype from the method and the processing
    """
    try:
        if (processing in PROCESSINGS_BIOTYPE):
            transform_biotype = PROCESSINGS_BIOTYPE[processing]
            if (re.search('->', transform_biotype)):
                modes = transform_biotype.split("->")
                biotype = re.sub(modes[0], modes[1], TYPES_BIOTYPE[METHODS_TYPE[method.lower()]])
            else:
                biotype = PROCESSINGS_BIOTYPE[processing]
        elif (method == "uniform"):
            biotype = "Discrete Copy number data" # TODO Continuous is better for grouping but posses problems with glyphs
        else:
            biotype = TYPES_BIOTYPE[METHODS_TYPE[method.lower()]]
    except:
        biotype = TYPES_BIOTYPE["unknown"]
        warn("Biotype of " + method.lower() + " is unknown")
    return biotype

def getDataType(method, processing="raw"):
    """
        Get the type of data from the method and the processing
    """
    try:
        datatype = METHODS_TYPE[method.lower()]
    except:
        datatype = "unknown"
    return datatype
