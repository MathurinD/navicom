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

DEBUG_NAVIDATA = True

# Constants related to NaviCell
MAX_GLYPHS = 5
GLYPH_TYPES = ["color", "size", "shape"]
# Identify the categories of cBioportal data as (aliases, biotype)
TYPES_SPEC = dict()
TYPES_SPEC["mRNA"] = (["mrna", "zscores", "mrna_median_zscores", "rna_seq_mrna_median_zscores", "rna_seq_mrna", "rna_seq_v2_mrna", "rna_seq_v2_mrna_median_zscores", "mrna_U133", "mrna_U133_zscores", "mrna_median", "mrna_zbynorm", "mrna_outliers", "rna_seq_rna", "mrna_znormal", "mrna_outlier", "mrna_merged_median_zscores"], "mRNA expression data")
TYPES_SPEC["dCNA"] = (["gistic", "cna", "cna_rae", "cna_consensus", "snp-fasst2"], "Discrete Copy number data")
TYPES_SPEC["cCNA"] = (["log2cna"], "Continuous copy number data")
TYPES_SPEC["methylation"] = (["methylation", "methylation_hm27", "methylation_hm450"], "mRNA expression data")
TYPES_SPEC["protein"] = (["RPPA_protein_level"], "protein level")
TYPES_SPEC["miRNA"] = (["mirna", "mirna_median_zscores"], "miRNA expression data")
TYPES_SPEC["mutations"] = (["mutations"], "Mutations")
TYPES_SPEC["unknown"] = (["unknown"], "mRNA expression data") # If the type of data cannot be identified, consider continuous data by default
METHODS_TYPE = dict()
TYPES_BIOTYPE = dict()
for bt in TYPES_SPEC:
    TYPES_BIOTYPE[bt] = TYPES_SPEC[bt][1]
    for cat in TYPES_SPEC[bt][0]:
        METHODS_TYPE[cat] = bt

# Inform what the processing does to the biotype, if -> it changes only some biotypes, if "something" it turns everything to something, see exportData
PROCESSINGS = ["raw", "moduleAverage", "pcaComp", "geoSmooth", "distribution", "colors"]
PROCESSINGS_BIOTYPE = {"moduleAverage":"Discrete->Continous", "pcaComp":"Color", "geoSmooth":"Discrete->Continuous"} 


class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name
    """
    def __init__(self, data, rows_list, columns_list, processing="raw", method="unknown", dType="data"):
        assert(len(data) == len(rows_list))
        for line in data:
            assert len(line) == len(columns_list), "Incorrect length of line : " + str(line) + ", length = " + str(len(line)) + ", expected " + str(len(columns_list))
        # Informations on the data TODO
        self.processing = processing
        self.method = method
        self.biotype = getBiotype(method, processing)
        # Initialise data and indexes
        self.data = np.array(data)
        self.rows = listToDictKeys(rows_list)
        self.rows_names = list(self.rows.keys())
        self.columns = listToDictKeys(columns_list)
        self.columns_names = list(self.columns.keys())
        if (dType == "annotations"):
            self.inColumns = "annotations"
            self.annotations = self.columns
            self.annotations_names = self.columns_names
            self.inRows = "samples"
            self.samples = self.rows
            self.samples_names = self.rows_names
        elif (dType == "data"):
            self.inRows = "genes"
            self.genes = self.rows
            self.genes_names = self.rows_names
            self.inColumns = "samples"
            self.samples = self.columns
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

    def __getitem__(self, index):
        if (isinstance(index, int)):
            return(self.data[index])
        elif (isinstance(index, str)):
            if (index in self.rows):
                return( NaviSlice(self.data[self.rows[index],:], self.columns) )
            elif (index in self.columns):
                return( NaviSlice(self.data[:,self.columns[index]], self.rows) )
            else:
                raise IndexError("'" + index + "' is neither a gene or sample name")
        elif (isinstance(index, list) or isinstance(index, tuple)):
            result = list()
            for ii in index:
                result.append(self[ii].data)
            result = np.array(result)
            if (index[0] in self.rows):
                assert np.all([idx in self.rows for idx in index]), "Not all index are gene names"
                genes = self.rows.copy()
                for gene in self.rows:
                    if (not gene in index):
                        genes.pop(gene)
                return( NaviData(result, genes, self.columns) )
            elif (index[0] in self.columns):
                assert np.all([idx in self.columns for idx in index]), "Not all index are sample names"
                samples = self.columns.copy()
                for sample in self.columns:
                    if(not sample in index):
                        samples.pop(sample)
                return( NaviData(result.transpose(), self.rows, samples) )

    def __iter__(self, by="genes"):
        if (not by in ["genes", "samples"]):
            raise ValueError("'by' must be in ['genes', 'samples']")
        self.iter_mode = by
        self.index = 0
        return(self)

    def __next__(self):
        try:
            if (self.iter_mode == "genes"):
                key = list(self.rows.keys())[self.index]
                result = [key] + list(self.data[self.index,:])
            elif (self.iter_mode == "samples"):
                key = list(self.columns.keys())[self.index]
                result = [key] + list(self.data[:,self.index])
        except IndexError:
            raise StopIteration
        self.index += 1
        return(result)

    def __repr__(self):
        rpr = "NaviData array with " + str(len(self.rows)) + " " + self.inRows + " and "
        rpr += str(len(self.columns)) + " " + self.inColumns
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
            ret = "GENE\t" + '\t'.join(self.columns_names) + "\n"
            for line in self:
                if (line[0] in hugo_map or hugo_map == ""):
                    ret += buildLine(line)
        else:
            raise ValueError("Data type must be 'data' or 'annotations'")

        return("@DATA\n" + ret)

    def exportToNaviCell(self, nv, biotype, dataName):
        """
        Export data to a NaviCell map
        """
        nv.importDatatables(self.makeData(nv.getHugoList()), dataName, biotype) # Remove name once in self
        #nv.importDatatables(self.makeData(nv.getHugoList()), dataName, self.biotype)

    def saveData(self, baseName="", mode="w"):
        """
        Save the NaviData datas in a file that can be used in NaviCell, but can also be loaded as NaviData
        """
        fname = baseName + "[" + self.processing + "]" + "[" + self.method + "].ncd"
        if (mode == "a"):
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
        self.samplesPerCategory = dict()
        modified_data = list()
        modified_annot = list()
        for annot in self.annotations:
            self.categoriesPerAnnotation[annot] = list()
            self.samplesPerCategory[annot] = dict()
            reduced = False
            # Reduce the annotations set if they are continuous integers with too many values
            if (len(np.unique(self[annot].data)) > MAX_GROUPS):
                try:
                    values = [float(value) for value in np.unique(self[annot].data)]
                    # Define the new annotations
                    min_value = min(values)
                    max_value = max(values)
                    step = (max_value-min_value)/MAX_GROUPS
                    new_values = [signif(cat) for cat in np.arange(min_value, max_value+step/2, 1.01*step)]
                    new_categories = [str(new_values[icat])+"<X<"+str(new_values[icat+1]) for icat in range(len(new_values)-1)] + ["NaN"]
                    self.categoriesPerAnnotation[annot] = new_categories
                    # Attribute the new annotations to samples
                    modified_data.append(self[annot].data.copy())
                    modified_annot.append(annot)
                    for cat in new_categories:
                        self.samplesPerCategory[annot][cat] = list()
                    for sample in self.samples:
                        annot_value = float(self[annot][sample])
                        if (np.isnan(annot_value)):
                            self.samplesPerCategory[annot]["NaN"].append(sample)
                        else:
                            icat = 0
                            while (signif(annot_value) >= new_values[icat+1]):
                                icat += 1
                            self.samplesPerCategory[annot][new_categories[icat]].append(sample)
                            self[annot][sample] = new_categories[icat]
                    reduced = True
                except ValueError: 
                    reduced = False
            # Non numeric values or with few enough levels are not modified and simply indexed
            if (not reduced):
                for sample in self.samples:
                    annot_value = self[annot][sample]
                    if (not annot_value in self.categoriesPerAnnotation[annot]):
                        self.categoriesPerAnnotation[annot].append(annot_value)
                        self.samplesPerCategory[annot][annot_value] = [sample]
                    else:
                        self.samplesPerCategory[annot][annot_value].append(sample)

        if (DEBUG_NAVIDATA):
            print("Modified annotations:" + str(modified_annot))
        if (len(modified_annot) > 0):
            self.old_annots = NaviData(modified_data, modified_annot, self.rows, "old_annotations")
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
    return(round(x, -int(math.log10(x))+(n-1) ))

def getBiotype(method, processing="raw"):
    """
    Get the biotype from the method and the processing
    """
    if (processing in PROCESSINGS_BIOTYPE):
        transform_biotype = PROCESSINGS_BIOTYPE[processing]
        if (re.search('->', transform_biotype)):
            modes = transform_biotype.split("->")
            biotype = re.sub(modes[0], modes[1], TYPES_BIOTYPE[METHODS_TYPE[method.lower()]])
        else:
            biotype = PROCESSINGS_BIOTYPE[processing]
    else:
        biotype = TYPES_BIOTYPE[METHODS_TYPE[method.lower()]]
    return biotype

