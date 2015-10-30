"""
#\mainpage navicom package documentation
The navicom package intends to provide standard methods to visualise high-throughput data in NaviCell web service. It also provide several processing method to get some extra meaning out of the data.

\author Mathurin Dorel
\par License:\n Lesser GNU Public License.

\section start Getting started

The communication with NaviCell web service, and data processing are performed by the \ref navicom::navicom::NaviCom "NaviCom" class, which can be initialised with a simple NaviCell map URL.
\code
nc = NaviCom(map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php')
\endcode
Data and annotations can then be loaded in the \ref navicom::navicom::NaviCom "NaviCom" object.
\code
nc.loadData("data/Ovarian_Serous_Cystadenocarcinoma_TCGA_Nature_2011.txt")
\endcode
Two formats are accepted:
    \li a matrix format where data are represented as matrix, with explicit rows and columns names. The first row has to start by GENE (for data) or NAME (for annotations);
    \li a set of matrix format where each matrix is separated by a header line expliciting the type of data and starting with M (data) or ANNOTATIONS (annotations).

After being loaded, data and annotations can be easily controlled:
\code
nc.listData()
nc.listAnnotations()
\endcode

\subsection display_config Display configuration

The displays in the NaviCell map can be easily configured with the \ref navicom::displayConfig::DisplayConfig "DisplayConfig" class. It can be provided to the \ref navicom::navicom::NaviCom "NaviCom" class.
\code
display_config = DisplayConfig(5, zero_color="000000", na_color="ffffff", na_size=0)
nc = NaviCom(map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php', display_config)
\endcode

\subsection extra_data Adding extra data

Data in the \ref navicom::navicom::NaviCom "NaviCom" class are represented in the \ref navicom::navidata::NaviData "NaviData" format, which is a matrix of data that can be indexed by row and column names. \ref navicom::navidata::NaviData "NaviData" objects can be created independently and integrated to a \ref navicom::navicom::NaviCom "NaviCom" object:
\code
extra_data = NaviData(data_matrix, row_names, col_names, method, processing)
nc.bindNaviData(extra_data, "extra_data")
\endcode
Each datatable is identified by its method, which is the biological or computationnal method used to generate the data, and its processing, which is a NaviCell map related processing to tweak the visualisation of the data. The default processing for unprocessed data is 'raw".

\section visualisation Data visualisation

The \ref navicom::navicom::NaviCom "NaviCom" class provides several method to display the data in NaviCell:
    \li \b display Generic display function to perform any kind of personnalised display
    \li \b displayOmics Display -omics data as map staining with extra information or data on using the other display modes
    \li \b completeDisplay Display all types of data available using map staining for expression, heatmap for copy number and glyph size for the others
    \li \b displayMutations Display mutations data as glyphs
    \li \b displayMutationsWithGenomics Display expression data as map staining, copy number as heatmaps and mutations as glyphs
    \li \b displayExpression Display expession data as map staining
    \li \b displayExpressionWithMutations Display expession data as map staining and mutations as glyphs
    \li \b displayExpressionWithCopyNumber Display expession data as map staining and copy number as heatmap
    \li \b displayExpressionWithProteomics Display expession data as map staining and proteomics levels as barplot
    \li \b displayExpressionWithMethylation Display expession data as map staining and methylation levels as barplot
    \li \b displayExpressionWithmiRNA Display expession data as map staining and miRNA levels as barplot

"""
