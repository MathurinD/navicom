# Import data from a file and create analyzed sets to export
from curie.navicell import *
import numpy as np
import re
from time import sleep
from collections import OrderedDict as oDict
from pprint import pprint

DEBUG_NAVICOM = True
VERBOSE_NAVICOM = True

# Constants related to NaviCell
MAX_GLYPHS = 5
GLYPH_TYPES = ["color", "size", "shape"]
# Identify the categories of cBioportal data as (aliases, biotype)
category_biotypes = dict()
category_biotypes["mRNA"] = (["mrna", "zscores", "mrna_median_Zscores", "rna_seq_mrna_median_Zscores", "rna_seq_mrna", "rna_seq_v2_mrna", "rna_seq_v2_mrna_median_Zscores", "mrna_U133", "mrna_U133_Zscores", "mrna_median", "mrna_zbynorm", "mrna_outliers", "rna_seq_rna", "mrna_znormal", "mrna_outlier"], "mRNA expression data")
category_biotypes["dCNA"] = (["gistic", "cna", "CNA", "cna_rae", "cna_consensus", "SNP-FASST2"], "Discrete Copy number data")
category_biotypes["cCNA"] = (["log2CNA"], "Continuous copy number data")
category_biotypes["METHYLATION"] = (["methylation", "methylation_hm27", "methylation_hm450"], "Methylation data")
category_biotypes["PROT"] = (["RPPA_protein_level"], "protein level")
category_biotypes["miRNA"] = (["mirna", "mirna_median_Zscores", "mrna_merged_median_Zscores"], "miRNA expression data")
category_biotypes["mutations"] = (["mutations"], "Mutations")
category = dict()
biotypes = dict()
for bt in category_biotypes:
    biotypes[bt] = category_biotypes[bt][1]
    for cat in category_biotypes[bt][0]:
        category[cat] = bt

PROCESSINGS = ["raw", "moduleAverage", "pcaComp", "geoSmooth"]
PROCESSINGS_BIOTYPES = {"moduleAverage"="Discrete->Continous", "pcaComp"="Color", "geoSmooth"="Discrete->Continuous"} # Inform what the processing does to the biotype, if -> it changes only some biotypes, if "something" it turns everything to something

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
                print("\t" + processing + " " + dname + ": " + biotypes[category[dname]])

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
                self.annotations = NaviData(annot["annot"], annot["samples"], annot["names"], dType="annotations")

    def defineUniformData(self, samples, genes):
        self.data["uniform"] = NaviData( np.array([[1] * len(samples) for nn in genes]), genes, samples )
        self.exportData("uniform")

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

    def display(self, perform_list, samples="all: 1.0", colors="", module=''):
        """
        Display data on the NaviCell map
        Args :
            perform_list (list of 2-tuples): each tuple must contain the name of the data to display and the mode of display ("glyphN_(color|size|shape)", "barplot", "heatmap" or "map_staining"). Barplots and heatmaps cannot be displayed simultaneously. Several data types can be specified for heatmaps. Specifying "glyph" (without number) will automatically select a new glyph for each data using the same properties (shape, color or size) in glyphs (maximum of 5 glyphs).
            colors : range of colors to use (NOT IMPLEMENTED YET)
            samples (str or list of str) : Samples to use. Only the first sample is used for glyphs and map staining, all samples from the list are used for heatmaps and barplots. Use 'all_samples' to use all samples or ['annot1:...:annotn', 'all_groups'] to use all groups corresponding to the combinations of annot1...annotn.
        """
        assert isinstance(perform_list, list), "perform list must be a list"
        assert isinstance(perform_list[0], tuple) and len(perform_list[0]) == 2, "perform list must be a list of 2-tuples"
        self.checkBrowser()
        self.exportAnnotations()
        self.resetDisplay()

        # Selection of samples
        all_samples = False
        all_groups = False
        if (isinstance(samples, str)):
            one_sample = samples
            samples = [samples]
        elif (isinstance(samples, list)):
            one_sample = samples[0]
        # Select the groups that must be selected to produce the composite groups required and control that all groups are compatible (because lower order composition are not generated)
        if (samples[0] == "all" or samples[0] == "all_samples" or samples[0] == "samples"):
            all_samples = True
        elif (len(samples) > 1 and (samples[1] == "all_groups" or samples[1] == "groups")):
            all_groups = True
        rGroups = 0 
        groups_list = []
        first_groups = True
        self.resetAnnotations()
        for sample in samples:
            nGroups = 0
            selections = sample.split(";")
            groups = []
            for select in selections:
                groups.append(select.split(":")[0].strip())
            # TODO Control that groups exist
            # TODO Check how to use combined groups
            if (len(groups) > 1 or groups[0] in self.annotations.annotations): # Skip the loop for samples
                if (first_groups):
                    first_groups = False
                    for group in groups:
                        if (group in self.annotations.annotations):
                            if (DEBUG_NAVICOM):
                                print("Selecting " + group)
                            #self.selectAnnotations(group) # Useless as long as annotations are removed on data import
                            groups_list.append(group)
                            nGroups += 1
                            rGroups += 1
                else:
                    for group in groups:
                        if (group in self.annotations.annotations):
                            assert group in groups_list, "Groups combinations are not compatibles"
                            nGroups += 1
                if (nGroups != rGroups and rGroups != 0 and nGroups != 0):
                    raise ValueError("Groups combinations are not compatible")

        # Control that the user does not try to display to many data or use several time the same display
        glyph = {gtype:[False] * MAX_GLYPHS for gtype in GLYPH_TYPES}
        if (len(samples) == 1):
            glyph_samples = samples * MAX_GLYPHS
            sample_for_glyph = [True] * MAX_GLYPHS
        else:
            glyph_samples = samples
            sample_for_glyph = [True] * len(samples)
            while (len(glyph_samples) < MAX_GLYPHS):
                glyph_samples.append(None)
                sample_for_glyph.append(False)
        glyph_set = False
        heatmap = False
        barplot = False
        barplot_data = ""
        map_staining = False

        # Preprocess the perform list to get valid data_name, and export data that have not been exported yet
        for perf_id in range(len(perform_list)):
            data_name = perform_list[perf_id][0]
            if (isinstance(data_name, str)):
                if (not re.match("_", data_name)):
                    data_name = data_name + "_raw"
                data_name = self.associated_data[data_name]
            processing = data_name[0]
            method = data_name[1]
            assert processing in self.processings, "Processing " + processing + " does not exist"
            self.exportData(method, processing)
            perform_list[perf_id] = (self.data_names[processing][method], perform_list[perf_id][1])

        self.selectAnnotations(groups_list) # Write annotations AFTER the export
        # Perform the display depending of the selected mode
        for data_name, dmode in perform_list:
            dmode = dmode.lower()
            if (re.match("^(glyph|color|size|shape)", dmode)):
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
                
                if (not glyph_number in range(1, MAX_GLYPHS+1)):
                    raise ValueError("Glyph number must be in [1," + str(MAX_GLYPHS) + "]")
                if (not glyph_type in GLYPH_TYPES):
                    raise ValueError("Glyph type must be one of " + str(GLYPH_TYPES))
                if (glyph[glyph_type][glyph_number-1]):
                    raise ValueError(glyph_type + " for glyph " + str(glyph_number) + " has already been specified")
                glyph[glyph_type][glyph_number-1] = True
                if (not glyph_samples[glyph_number-1]):
                    raise ValueError("Incorrect glyph number : " + str(glyph_number) + ", only " + str(len(samples)) + " have been given")

                cmd="self.nv.glyphEditorSelect" + glyph_type.capitalize() + "Datatable(" + module +  ", " + str(glyph_number) + ", '" + data_name + "')"
                print(cmd)
                exec(cmd)
            elif (re.match("map_?staining", dmode)):
                if (not map_staining):
                    self.nv.mapStainingEditorSelectDatatable(module, data_name)
                    self.nv.mapStainingEditorSelectSample(module, one_sample)
                    self.nv.mapStainingEditorApply(module)
                    map_staining = True
                else:
                    raise ValueError("Map staining can only be applied once, use a separate call to the display function to change map staining")
            elif (re.match("heatmap", dmode)):
                if (barplot):
                    raise ValueError("Heatmaps and barplots cannot be applied simultaneously, use a separate call to the display function to perform the heatmap")
                else:
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
            elif (re.match("barplot", dmode)):
                if (heatmap):
                    raise ValueError("Heatmaps and barplots cannot be applied simultaneously, use a separate call to the display function to perform the barplot")
                else:
                    # Check that it does not try to add new data, and simply adds samples
                    if (barplot and data_name != barplot_data and not re.search("all|groups|samples", data_name)):
                        raise ValueError("Barplot has already been set with different data, use a separate call to the display function to perform another barplot")
                    else:
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
            else:
                raise ValueError(dmode + " drawing mode does not exist")

        # Check that datatables are selected for all glyphs features (until default has been added)
        # or complete, then apply the glyphs configuration
        if (glyph_set):
            for glyph_id in range(MAX_GLYPHS):
                nsets = sum(1 for cs in GLYPH_TYPES if glyph[cs][glyph_id])
                if (nsets > 0):
                    if (not sample_for_glyph[glyph_id]):
                        raise ValueError("No sample has been attributed to glyph " + str(glyph_id+1))
                    else:
                        self.nv.glyphEditorSelectSample(module, glyph_id+1, glyph_samples[glyph_id])
                    if (not glyph["color"][glyph_id]):
                        self.nv.glyphEditorSelectColorDatatable(module, glyph_id+1, "uniform")
                    if (not glyph["shape"][glyph_id]):
                       self.nv.glyphEditorSelectShapeDatatable(module, glyph_id+1, "uniform")
                    if (not glyph["size"][glyph_id]):
                        self.nv.glyphEditorSelectSizeDatatable(module, glyph_id+1, "uniform")
                    self.nv.glyphEditorApply(module, glyph_id+1)
        if (barplot):
            self.nv.barplotEditorApply(module)
        if (heatmap):
            self.nv.heatmapEditorApply(module)

    #def displayGroups(self, groups, combine=T, method): # method in ["barplot", "heatmap", "glyph_TYPE"]
    
    #def addDisplay(self, perform_list, samples="all", colors=""):

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
            if (method in self.data[processing] and method in self.methodBiotype):
                if (not self.exported_data[processing][method]):
                    name = self.nameData(method, processing, name)
                    # Processing turn discrete data into continuous or color data
                    if (processing in self.methodBiotype):
                        biotype = self.methodBiotype[processing]
                    else:
                        biotype = self.methodBiotype[method]
                    self.nv.importDatatables(self.data[processing][method].makeData(self.nv.getHugoList()), name, biotype)
                    self.exported_data[processing][method] = True
                    done_export = True
                elif (VERBOSE_NAVICOM):
                    print(method + " data with " + processing + " processing has already been exported")
            elif (method == "uniform"): # Uniform data for glyphs
                self.nv.importDatatables(self.data["uniform"].makeData(self.nv.getHugoList()), "uniform", "Discrete Copy number data") # Continuous is better for grouping but posses problems with glyphs
                name = "uniform"
                done_export = True
            else:
                raise KeyError("Method " + method + " with processing " + processing + " does not exist")
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

    def displayMethylome(self, background="mRNA", processing="raw", patient="all: 1.0"):
        """
            Display the methylation data as glyphs on the NaviCell map, with mRNA expression of gene CNV as map staining
            Args:
                background (str) : should genes, mRNA or no data be used for the map staining
                processing (str) : should the processed data be used
        """
        mrna_alias = ["MRNA"] # TODO Define what to put here
        gene_alias = ["CNV", "CNA", ""]
        assert background.upper in mrna_alias + gene_alias + ["NO", ""], "Select either genes, mRNA or no data for the map staining"
        # Display methylation as glyphs
        for method in self.data[processing]:
            if (re.match("methylation", method.lower())):
                self.exportData(self, method, processing)
        # Display mRNA or gene data as map staining
        data = self.data[processing][background]



class NaviData():
    """
    Custom class to store the data and be able to access rows and columns by name
    """

    def __init__(self, data, genes_list, samples_list, dType="data"):
        assert(len(data) == len(genes_list))
        for line in data:
            assert len(line) == len(samples_list), "Incorrect length of line : " + line + ", length = " + str(len(line))
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
    return(re.sub("NaN|nan", "NA", ret[:-1]) + "\n")




