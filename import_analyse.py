# Import data from a file and create analyzed sets to export
from curie.navicell import *
import numpy as np
import re
from time import sleep
from collections import OrderedDict as oDict
from pprint import pprint
import math

DEBUG_NAVICOM = True
VERBOSE_NAVICOM = True

# Constants related to NaviCell
MAX_GLYPHS = 5
GLYPH_TYPES = ["color", "size", "shape"]
# Identify the categories of cBioportal data as (aliases, biotype)
TYPES_SPEC = dict()
TYPES_SPEC["mRNA"] = (["mrna", "zscores", "mrna_median_Zscores", "rna_seq_mrna_median_Zscores", "rna_seq_mrna", "rna_seq_v2_mrna", "rna_seq_v2_mrna_median_Zscores", "mrna_U133", "mrna_U133_Zscores", "mrna_median", "mrna_zbynorm", "mrna_outliers", "rna_seq_rna", "mrna_znormal", "mrna_outlier"], "mRNA expression data")
TYPES_SPEC["dCNA"] = (["gistic", "cna", "CNA", "cna_rae", "cna_consensus", "SNP-FASST2"], "Discrete Copy number data")
TYPES_SPEC["cCNA"] = (["log2CNA"], "Continuous copy number data")
TYPES_SPEC["methylation"] = (["methylation", "methylation_hm27", "methylation_hm450"], "mRNA expression data")
TYPES_SPEC["PROT"] = (["RPPA_protein_level"], "protein level")
TYPES_SPEC["miRNA"] = (["mirna", "mirna_median_Zscores", "mrna_merged_median_Zscores"], "miRNA expression data")
TYPES_SPEC["mutations"] = (["mutations"], "Mutations")
METHODS_TYPE = dict()
TYPES_BIOTYPE = dict()
for bt in TYPES_SPEC:
    TYPES_BIOTYPE[bt] = TYPES_SPEC[bt][1]
    for cat in TYPES_SPEC[bt][0]:
        METHODS_TYPE[cat] = bt

# Inform what the processing does to the biotype, if -> it changes only some biotypes, if "something" it turns everything to something, see exportData
PROCESSINGS = ["raw", "moduleAverage", "pcaComp", "geoSmooth"]
PROCESSINGS_BIOTYPE = {"moduleAverage":"Discrete->Continous", "pcaComp":"Color", "geoSmooth":"Discrete->Continuous"} 

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
    """
    NaviComm class to handle data and display them in a standardized way on NaviCell maps
    """

    def __init__(self, fname="data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt", map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php', modules_dict=""):
        # Build options for the navicell connexion
        options = Options()
        options.map_url = map_url
        idx = options.map_url.find('/navicell/')
        options.proxy_url = options.map_url[0:idx] + '/cgi-bin/nv_proxy.php'
        options.browser_command = "firefox %s" # TODO Add user control
        self.nv = NaviCell(options)
        # NaviCell export control
        self.exported_annotations = False
        self.browser_opened = False
        self.biotypes = dict()
        self.methodBiotype = {"gistic":"Discrete Copy number data", "log2CNA":"Continuous copy number data", "rna_seq_mrna":"mRNA expression data", "moduleAverage":"Continuous copy number data"} # TODO Delete and change to biotypes[category]
        # Data, indexed by processing then type of data
        self.processings = ["raw", "moduleAverage", "pcaComp", "geoSmooth"]
        self.data = dict() # All data
        self.exported_data = dict() # Whether data have been exported yet
        self.data_names = dict() # Name of the exported data
        self.associated_data = dict() # Processing and method associated to each name
        for processing in self.processings:
            self.data[processing] = dict()
            self.exported_data[processing] = dict()
            self.data_names[processing] = dict()
        self.exported_data["uniform"] = False
        self.annotations = dict() # Annotations of the samples
        self.modules = dict() # Composition of each module
        self.associated_modules = dict() # Number of modules each gene belong to
        if (fname != ""):
            self.loadData(fname)
            self.defineModules(modules_dict)
        # Remember how many samples and datatables were selected in NaviCell in the heatmap and barplot editors
        self.hsid = 0
        self.hdid = 0
        self.bid = 0

    def dataAvailable(self):
        print("Data available :")
        for processing in self.processings:
            for dname in self.data[processing]:
                if (dname in METHODS_TYPE):
                    print("\t" + processing + " " + dname + ": " + METHODS_TYPE[dname] + " (biotype:" + TYPES_BIOTYPE[METHODS_TYPE[dname]] + ")")
                else:
                    print("\t" + processing + " " + dname)

    def __repr__(self):
        rpr = "NaviCom object with " + str(len(self.data)) + " types of data:\n"
        for method in self.data:
            rpr += method + ": " + self.data[method] + "\n"
        rpr += "and " + str(len(self.moduleAverage)) + " modules average:\n"
        for method in self.moduleAverage:
            rpr += method + ": " + self.moduleAverage[method] + "\n"
        return(repr)
    
    def nameData(self, method, processing="raw", name=""):
        if (name == ""):
            name = method + "_" + processing
        self.data_names[processing][method] = name
        self.associated_data[name] = (processing, method)
        return (name)

    def getDataName(self, data_name):
        if (isinstance(data_name, str))
            if (data_name in self.associated_data):
                return(data_name)
            elif (data_name + "_raw" in self.associated_data):
                return(data_name + "_raw")
        elif (isinstance(data_name, tuple) and len(data_name == 2)):
            if (data_name[0] in self.processings):
                return(self.data_name[data_name[0]][data_name[1]])
            elif (data_name[1] in self.processings):
                return(self.data_name[data_name[1]][data_name[0]])
        raise ValueError("Invalid name for data: " + str(data_name))

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
                self.nameData(method, "raw")
                if (not "uniform" in self.data):
                    self.defineUniformData(profile_data["samples"], profile_data["genes"])
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
                self.annotations = NaviAnnotations(annot["annot"], annot["samples"], annot["names"], dType="annotations")

    def defineUniformData(self, samples, genes):
        self.data["uniform"] = NaviData( np.array([[1] * len(samples) for nn in genes]), genes, samples )
        self.exportData("uniform")

    def newProcessedData(self, processing, method, data):
        assert processing in self.processings, "Processing " + processing + " is not handled"
        self.data[processing][method] = data
        self.exported_data[processing][method] = False
        self.nameData(method, processing)

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
                            self.modules[module_name] += ll[1:]
                        else:
                            self.modules[module_name] += ll[2:]
                        # Count the number of modules each gene belong to
                        for gene in self.modules[module_name]:
                            try:
                                self.associated_modules[gene].append(module_name)
                            except KeyError:
                                self.associated_modules[gene] = [module_name]
        """
        # Only keep genes with data
        for module in self.modules:
            keep = list()
            for gene in self.modules[module]:
                if (not gene in self.genes_list):
                    self.modules[module].remove(gene)
        """

    def averageModule(self, method):
        """
        Perform module averaging for every modules for one data type
        """
        assert method in self.data["raw"], "This type of data is not present"
        assert len(self.modules)>0, "No module have been defined"

        # Calculate average expression for each module
        data = self.data["raw"][method]
        samples = list(data.samples.keys())
        module_expression = dict()
        for module in self.modules:
            module_expression[module] = [0 for sample in samples]
            non_nan = np.array([0 for sample in samples])
            no_data = list()
            for gene in self.modules[module]:
                try:
                    not_nan = np.array([int(not np.isnan(dd)) for dd in data[gene].data])
                    non_nan += not_nan
                    non_nan_data = data[gene].data * not_nan
                    non_nan_data[np.isnan(non_nan_data)] = 0
                    module_expression[module] += non_nan_data
                except IndexError:
                    no_data += [gene]
            for gene in no_data:
                #self.modules[module].remove(gene)
                if (VERBOSE_NAVICOM):
                    print(gene + " from module " + module + " has no " + method + " data")
            module_expression[module] /= non_nan

        # Calculate average module expression for each gene
        gene_module_average = list()
        for gene in data.genes:
            if gene in self.associated_modules:
                gene_module_average.append(np.array([0. for sample in data.samples]))
                for module in self.associated_modules[gene]:
                    gene_module_average[-1] += module_expression[module] / len(self.associated_modules[gene])
            else:
                gene_module_average.append(np.array(list(data[gene])))
                print(gene)

        # Put the averaging in a NaviData structure
        #self.data["moduleAverage"][method] = NaviData(list(module_expression.values()), list(self.modules.keys()), samples) # Usefull if NaviCell allow modules values one day
        self.newProcessedData("moduleAverage", method, NaviData(gene_module_average, list(data.genes), samples))

    def pcaComp(self, method, colors=["red", "green", "blue"]):
        """
        Run pca on the data and create a color matrix with the 3 principal components in the three main colors
        """
        print("Not implemented yet")

    def exportData(self, method, processing="raw", name=""):
        """
        Export data to NaviCell, can be processed data

        Args:
            method (str) : name of the method to export
            processing (str) : "" to export raw data, processing method to export processed data. See 'averageModule' and 'pcaComponent'
        """
        self.checkBrowser() # TODO Perform processing if necessary
        done_export = False

        if (processing in self.processings):
            if (method in self.data[processing]):
                if(method in METHODS_TYPE):
                    if (not self.exported_data[processing][method]):
                        name = self.nameData(method, processing, name)
                        # Processing change the type of data, like discrete data into continuous, or anything to color data
                        if (processing in PROCESSINGS_BIOTYPE):
                            transform_biotype = PROCESSINGS_BIOTYPE[processing]
                            if (re.search('->', transform_biotype)):
                                modes = transform_biotype.split("->")
                                biotype = re.sub(modes[0], modes[1], TYPES_BIOTYPE[processing])
                            else:
                                biotype = PROCESSINGS_BIOTYPE[processing]
                        else:
                            biotype = TYPES_BIOTYPE[METHODS_TYPE[method]]
                        self.nv.importDatatables(self.data[processing][method].makeData(self.nv.getHugoList()), name, biotype)
                        self.exported_data[processing][method] = True
                        done_export = True
                else:
                    raise ValueError("Biotype of '" + method + "' is unknown")
            elif (method == "uniform"): # Uniform data for glyphs
                if (not self.exported_data["uniform"]):
                    self.nv.importDatatables(self.data["uniform"].makeData(self.nv.getHugoList()), "uniform", "Discrete Copy number data") # Continuous is better for grouping but posses problems with glyphs
                    name = "uniform"
                    done_export = True
                    self.exported_data["uniform"] = True
            else:
                raise KeyError("Method '" + method + "' with processing '" + processing + "' does not exist")
        else:
            raise KeyError("Processing " + processing + " does not exist")

        # Sleep to avoid errors due to the fact that the loading by NaviCell is asynchronous
        # TODO Remove it when the python API receives signal
        if (done_export):
            print("Exporting " + name + " to NaviCell...")

    def checkBrowser(self):
        """
        Check if the browser is opened, and open it if it is not
        """
        if (not self.browser_opened):
            print("Launching browser...")
            self.nv.launchBrowser()
            self.browser_opened = True
        return(self.browser_opened)

    def exportAnnotations(self):
        """
        Export samples annotations to NaviCell
        """
        self.checkBrowser()

        if (not self.exported_annotations):
            self.nv.sampleAnnotationImport(self.annotations.makeData())
            self.exported_annotations = True

    def display(self, perform_list, default_samples="all: 1.0", colors="", module='', reset=True):
        """
        Display data on the NaviCell map
        Args :
            perform_list (list of 2-tuples): each tuple must contain the name of the data to display and the mode of display ("glyphN_(color|size|shape)", "barplot", "heatmap" or "map_staining"). Barplots and heatmaps cannot be displayed simultaneously. Several data types can be specified for heatmaps. Specifying "glyph" (without number) will automatically select a new glyph for each data using the same properties (shape, color or size) in glyphs (maximum of 5 glyphs).
            colors : range of colors to use (NOT IMPLEMENTED YET)
            default_samples (str or list of str) : Samples to use. Only the first sample is used for glyphs and map staining, all default_samples from the list are used for heatmaps and barplots. Use 'all_samples' to use all default_samples or ['annot1:...:annotn', 'all_groups'] to use all groups corresponding to the combinations of annot1...annotn.
        """
        assert isinstance(perform_list, list), "perform list must be a list"
        assert isinstance(perform_list[0], tuple) and (len(perform_list[0]) == 2 or len(perform_list[0]) == 3), "perform list must be a list of (2/3)-tuples"
        self.checkBrowser()
        self.exportAnnotations()
        if (reset):
            self.resetDisplay()

        # Preprocess the perform list to get valid data_name, and export data that have not been exported yet
        for perf_id in range(len(perform_list)):
            data_name = perform_list[perf_id][0]
            perform = perform_list[perf_id]
            if (isinstance(data_name, str)):
                if (not data_name in self.associated_data):
                    data_name = data_name + "_raw"
                data_name = self.associated_data[data_name]
            processing = data_name[0]
            method = data_name[1]
            assert processing in self.processings, "Processing " + processing + " does not exist"
            self.exportData(method, processing)
            if (len(perform) >= 3 and perform[2] != default_samples):
                perform_list[perf_id] = (self.data_names[processing][method], perform[1], perform[2])
            else:
                perform_list[perf_id] = (self.data_names[processing][method], perform[1], '')
        self.exportData("uniform")

        # Control that the user does not try to display to many data or use several times the same display
        if (True):
            glyph = {gtype:[False] * MAX_GLYPHS for gtype in GLYPH_TYPES}
            sample_for_glyph = [False] * MAX_GLYPHS
            glyph_samples = [""] * MAX_GLYPHS
            glyph_set = False
            heatmap = False
            barplot = False
            barplot_data = ""
            map_staining = False
            default_samples = self.processSamples(default_samples)
            lastWasDefault = True
            valid_default = (len(default_samples) == 1 and default_samples != "all_groups" and default_samples != "all_samples")
        # Perform the display depending of the selected mode
        for perform in perform_list:
            all_samples = False
            all_groups = False
            data_name = perform[0]
            dmode = perform[1]
            dmode = dmode.lower()
            # Check groups in NaviCell and get a valid list of samples, reload default if not the last used
            if (perform[2] == '' and not lastWasDefault):
                samples = self.processSamples(default_samples)
                lastWasDefault = True
            else:
                samples = self.processSamples(perform[2])
                lastWasDefault = False
            if (samples == "all_groups"):
                all_groups = True
            elif (samples == "all_samples"):
                all_samples = True

            if (re.search("^(glyph|color|size|shape)", dmode)):
                glyph_set = True
                # Extract the glyph id and the setup
                parse_setup = dmode.split("_")
                glyph_setup = parse_setup
                if (len(parse_setup) == 2):
                    try:
                        glyph_setup = [parse_setup[1] + str(int(parse_setup[0][-1]))]
                    except ValueError:
                        glyph_setup = [parse_setup[1]]
                elif (len(parse_setup) != 1):
                    raise ValueError("Glyph specification '" + dmode + "' incorrect")
                try:
                    glyph_number = int(glyph_setup[0][-1])
                    glyph_type = glyph_setup[0][:-1]
                except ValueError:
                    glyph_type = glyph_setup[0]
                    glyph_number = 1
                    while (glyph[glyph_type][glyph_number-1]):
                        glyph_number += 1
                glyph_id = glyph_number - 1
                
                if (not glyph_number in range(1, MAX_GLYPHS+1)):
                    raise ValueError("Glyph number must be in [1," + str(MAX_GLYPHS) + "]")
                if (not glyph_type in GLYPH_TYPES):
                    raise ValueError("Glyph type must be one of " + str(GLYPH_TYPES))
                if (glyph[glyph_type][glyph_number-1]):
                    raise ValueError(glyph_type + " for glyph " + str(glyph_number) + " has already been specified")
                glyph[glyph_type][glyph_number-1] = True
                if (len(samples) != 1 or samples == "all_groups" or samples == "all_samples"):
                    raise ValueError("Only one group or sample can be used for glyphs")

                if (lastWasDefault):
                    pass
                elif (glyph_samples[glyph_id] == ""):
                    glyph_samples[glyph_id] = samples[0]
                    sample_for_glyph[glyph_id] = True
                elif (glyph_samples[glyph_id] != samples[0]):
                    raise ValueError("Only one sample can be specified per glyph")

                cmd="self.nv.glyphEditorSelect" + glyph_type.capitalize() + "Datatable('" + module +  "', " + str(glyph_number) + ", '" + data_name + "')"
                if (DEBUG_NAVICOM):
                    print(cmd)
                exec(cmd)
            elif (re.search("map_?staining", dmode)):
                assert valid_default, "Only one sample can be used for map staining"
                if (not map_staining):
                    self.nv.mapStainingEditorSelectDatatable(module, data_name)
                    self.nv.mapStainingEditorSelectSample(module, samples[0])
                    self.nv.mapStainingEditorApply(module)
                    map_staining = True
                else:
                    raise ValueError("Map staining can only be applied once, use a separate call to the display function to change map staining")
            elif (re.search("heatmap", dmode)):
                if (barplot):
                    raise ValueError("Heatmaps and barplots cannot be applied simultaneously, use a separate call to the display function to perform the heatmap")
                heatmap = True
                # Select data
                self.nv.heatmapEditorSelectDatatable(module, self.hdid, data_name)
                self.hdid += 1
                # Select samples
                if (all_samples):
                    self.nv.heatmapEditorAllSamples(module)
                elif (all_groups):
                    self.nv.heatmapEditorAllGroups(module)
                elif (self.hsid == 0):
                    for spl in samples:
                        self.nv.heatmapEditorSelectSample(module, self.hsid, spl)
                        self.hsid += 1
                    self.nv.heatmapEditorApply(module)
            elif (re.search("barplot", dmode)):
                if (heatmap):
                    raise ValueError("Heatmaps and barplots cannot be applied simultaneously, use a separate call to the display function to perform the barplot")
                # Check that it does not try to add new data, and simply adds samples # TODO Remove as samples cannot be added (resetDisplay at the beginning)
                if (barplot and data_name != barplot_data and not re.search("all|groups|samples", data_name)):
                    raise ValueError("Barplot has already been set with different data, use a separate call to the display function to perform another barplot")
                barplot = True
                barplot_data = data_name
                self.nv.barplotEditorSelectDatatable(module, data_name)
                # Select samples
                if (all_samples):
                    self.nv.barplotEditorAllSamples(module)
                elif (all_groups):
                    self.nv.barplotEditorAllGroups(module)
                elif (self.bid == 0):
                    for spl in samples:
                        self.nv.barplotEditorSelectSample(module, self.bid, spl)
                        self.bid += 1
                    self.nv.barplotEditorApply(module)
            else:
                raise ValueError("'" + dmode + "' drawing mode does not exist")

        # Check that datatables are selected for all glyphs features (until default has been added)
        # or complete, then apply the glyphs configuration
        default_samples = self.processSamples(default_samples)
        if (glyph_set):
            for glyph_id in range(MAX_GLYPHS):
                nsets = sum(1 for cs in GLYPH_TYPES if glyph[cs][glyph_id])
                if (nsets > 0):
                    if (sample_for_glyph[glyph_id]):
                        self.nv.glyphEditorSelectSample(module, glyph_id+1, self.processSamples(glyph_samples[glyph_id]))
                    elif (valid_default):
                        print("Using default sample for glyph " + str(glyph_id+1))
                        self.nv.glyphEditorSelectSample(module, glyph_id+1, default_samples[0])
                    else:
                        raise ValueError("No samples specified for glyph " + str(glyph_id+1) + " and default_samples is invalid")
                    if (not glyph["color"][glyph_id]):
                        self.nv.glyphEditorSelectColorDatatable(module, glyph_id+1, "uniform")
                    if (not glyph["shape"][glyph_id]):
                       self.nv.glyphEditorSelectShapeDatatable(module, glyph_id+1, "uniform")
                    if (not glyph["size"][glyph_id]):
                        self.nv.glyphEditorSelectSizeDatatable(module, glyph_id+1, "uniform")
                    self.nv.glyphEditorApply(module, glyph_id+1)

    def resetDisplay(self):
        """
        Reset the data and samples selections in NaviCell
        """
        for ii in range(self.bid):
            self.nv.barplotEditorSelectSample('', ii, '')
        self.nv.barplotEditorSelectDatatable('', '')
        self.nv.drawingConfigSelectBarplot('', False)

        for ii in range(self.hsid):
            self.nv.heatmapEditorSelectSample('', ii, '')
        for ii in range(self.hdid):
            self.nv.heatmapEditorSelectDatatable('', ii, '')
        self.nv.drawingConfigSelectHeatmap('', False)

        for gid in range(1, MAX_GLYPHS):
            for gt in GLYPH_TYPES:
                exec("self.nv.glyphEditorSelect" + gt.capitalize() + "Datatable('', " + str(gid) + ", '')")
            self.nv.glyphEditorSelectSample('', gid, '')
            self.nv.drawingConfigSelectGlyph('', gid, False)
        
        self.nv.drawingConfigSelectMapStaining('', False)
        self.nv.drawingConfigApply('')

        # Reset the counters
        self.bid = 0
        self.hsid = 0
        self.hdid = 0

    def resetAnnotations(self, module=''):
        for annot in self.annotations.annotations:
            self.nv.sampleAnnotationSelectAnnotation(module, annot, False)
        self.nv.sampleAnnotationApply(module)

    def selectAnnotations(self, annotations, module=''):
        if (isinstance(annotations, str)):
            self.nv.sampleAnnotationSelectAnnotation(module, annotations, True)
        elif (isinstance(annotations, list)):
            for annot in annotations:
                self.nv.sampleAnnotationSelectAnnotation(module, annot, True)
        else:
            raise ValueError("'annotations' must be a string or a list")
        self.nv.sampleAnnotationApply(module)

    def processSamples(self, current_samples):
        """
        Process a list of samples or groups to a list of samples/groups names exportable to NaviCell or to "all_groups"/"all_samples" for heatmap and barplot, and select the correct groups in NaviCell
        """

        self.resetAnnotations()
        # Make sure current_samples is a list
        all_groups = False
        if (isinstance(current_samples, str)):
            current_samples = [current_samples]

        # It is possible to select all samples or all groups on heatmap and barplot
        if (current_samples[0] == "all_samples" or current_samples[0] == "samples"):
            return "all_samples"
        elif (current_samples[0] == "all" or current_samples[0] == "all: 1.0" or current_samples[0] == ""):
            self.selectAnnotations("all")
            return ["all: 1.0"]
        elif (len(current_samples) > 1 and (current_samples[1] == "all_groups" or current_samples[1] == "groups")):
            all_groups = True

        rGroups = 0 # Number of annotations to select
        groups_list = []
        first_groups = True
        # Select the groups that must be selected to produce the composite groups required
        for sample in current_samples:
            nGroups = 0
            groups = processGroups()[0]
            # Check that all groups are compatible in the annotations selected (because lower order composition are not generated). No check for individual samples
            if (len(groups) > 1 or groups[0] in self.annotations.annotations):
                if (first_groups): # Select the set of annotations for the first group
                    for group in groups:
                        if (group in self.annotations.annotations):
                            if (DEBUG_NAVICOM):
                                print("Selecting " + group)
                            self.selectAnnotations(group) 
                            groups_list.append(group)
                            nGroups += 1
                            rGroups += 1
                            first_groups = False
                        else:
                            raise ValueError("Annotation " + group + " does not exist")
                else: # Control that all other groups are compatible
                    for group in groups:
                        if (group in self.annotations.annotations):
                            assert group in groups_list, "Groups combinations are not compatibles as " + group + " is not in " + str(groups_list)
                            nGroups += 1
                if (nGroups != rGroups and rGroups != 0 and nGroups != 0):
                    raise ValueError("Groups combinations are not compatible, different number of groups")

        if (all_groups):
            return "all_groups"
        return current_samples

    def processGroups(self, groupName):
        """
        Process a group selection string and return the names of the individual groups to select and the corresponding values selected.
        """
        selections = groupName.split(";")
        groups = list()
        values = list()
        for select in selections:
            subName = select.split(":")
            groups.append(subName[0].split())
            try:
                values.append(float(subName[1].split()))
            except ValueError:
                values.append(subName[1].split())
        return((groups, values))

    def displayMethylome(self, samples="all: 1.0", processing="raw", background="mRNA", methylation="glyph"):
        """
            Display the methylation data as glyphs or heatmap on the NaviCell map, with mRNA expression of gene CNV as map staining
            Args:
                background (str) : should genes, mRNA or no data be used for the map staining
                processing (str) : should the processed data be used
        """
        mrna_alias = ["MRNA"] # TODO Define what to put here
        gene_alias = ["CNV", "CNA", "CCNA", "CCNV", "DNA"]
        background = background.upper()
        assert background in mrna_alias + gene_alias + ["NO", ""], "Select either genes, mRNA or no data for the map staining"
        assert methylation in ["glyph", "glyphs", "heatmap", "size", "glyph_size"], "Cannot use " + methylation + " to display methylation data"
        if (methylation != "heatmap"): methylation="size"
        assert processing in self.data, "Processing " + processing + " does not exist"

        # Select all methylation data and display as heatmap
        disp_selection = list()
        for method in self.data[processing]:
            included = method in METHODS_TYPE
            if (re.search("methylation", method.lower()) or (included and METHODS_TYPE[method] == "methylation")):
                disp_selection.append((self.data_names[processing][method], methylation)) # TODO Change to barplot when several datatables can be used
        # Display mRNA or gene data as map staining
        if (background in mrna_alias):
            for mrna in TYPES_SPEC["mRNA"][0]:
                if (mrna in self.data[processing]):
                    disp_selection.append((self.data_names[processing][mrna], "map_staining"))
                    break
        elif (background in gene_alias):
            for gene in TYPES_SPEC["cCNA"]+TYPES_SPEC["dCNA"]:
                if (gene in self.data[processing]):
                    disp_selection.append((self.data_names[processing][mrna], "map_staining"))
                    break
        if (DEBUG_NAVICOM):
            print(disp_selection)
            print(samples)
        self.display(disp_selection, samples)

    def displayTranscriptome(self, dataName, group="all: 1.0", samplesDisplay="", samples=list(), nbOfSamples=10):
        """
        Display one transcriptome data as map staining, and optionnaly samples from the group as heatmap
        Args:
            - dataName (str or tuple): name or identifier of the data.
            - group (str): Identifier of the group to display
            - samplesDisplay (str): Channel where the individual samples should be displayed (heatmap or barplot)
            - samples (list or str): list of samples to display, or a string specifying how such a list should be built ('quantiles' or 'random')
            - nbOfSamples (int): number of individual samples to display, ignored it samples is a list
        """
        allowedDisplays = ["", "heatmap", "barplot"]
        assert samplesDisplay in allowedDisplays, "samplesDisplay must one of " + str(allowedDisplays)
        dataName = getDataName(dataName)
        self.display([(dataName, "map_staining")], group)
        if (samplesDisplay != ""):
            if (isinstance(samples, list)):
                self.display([(dataName, samplesDisplay)], samples, reset=False)
            elif (isinstance(samples, str)):
                if (samples == "quantiles"):
                    pass
                elif (samples == "random"):
                    pass
                else:
                    self.display([(dataName, samplesDisplay)], [samples], reset=False)

    def selectQuantiles(self, dataName, group, numberOfQuantiles):
        pass

class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name
    """

    def __init__(self, data, rows_list, columns_list, dType="data"):
        assert(len(data) == len(rows_list))
        for line in data:
            assert len(line) == len(columns_list), "Incorrect length of line : " + str(line) + ", length = " + str(len(line))
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
            
MAX_GROUPS = 7
class NaviAnnotations(NaviData):
    """
    Enhance NaviData to contain annotations and associate annotations values with samples
    """

    def __init__(self, data, rows_list, columns_list, dType="annotations"):
        NaviData.__init__(self, data, rows_list, columns_list, dType)
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
                        self.samplesPerCategory[cat] = list()
                    for sample in self.samples:
                        annot_value = float(self[annot][sample])
                        if (np.isnan(annot_value)):
                            self.samplesPerCategory["NaN"].append(sample)
                        else:
                            icat = 0
                            while (signif(annot_value) >= new_values[icat+1]):
                                icat += 1
                            self.samplesPerCategory[new_categories[icat]].append(sample)
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

        if (DEBUG_NAVICOM):
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


