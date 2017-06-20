This file contains a description of the contents of this folder,  instructions for how to run Variant View, the variant data file format Variant View accepts, and license information.

(I) Package Contents

In addition to this READ_ME.txt file, this package includes the following files:

1) VariantView.html
2) VariantView.js
3) legend.png
4) Knowledge_Bases/uniprotIDFTREFSEQ.dat
5) Knowledge_Bases/refGene.csv
6) data_files/demo_data_set.txt

Files 1) and 2) contain the Variant View source code.  File 3) is the legend for the Variant View tool. File 4) is a trimmed flat file dump of the UniProt database for the human species. The original file was downloaded from http://www.uniprot.org/.  File 5) is a flat file dump from the RefSeq database containing gene transcript information such as exon coordinates, gene names, and gene ids; the original file was downloaded from http://www.ncbi.nlm.nih.gov/refseq/. File 6) is an example variant data set for demonstrating tool use and the proper format for a variant data file. The sample variant data set is provided by Dr. Linda Chang and Dr. Gerben Duns. All patient ids have been sanitized by replacing them with fictional names. The dataset is a subset of the use cases discussed in the paper, containing only the genes already well known to be implicated in acute myeloid leukemia (AML).

(II) Running Variant View

There are two options for running Variant View:

Option 1: Web Server

If you have administrative access to a web server, you can place the contents of this package in a directory on your server, and point your browser to the address associated with that directory plus "VariantView.html". For example, http://www.myexamplepage.com/VariantView.html

You can substitute your own dataset for the provided example by changing the contents of the demo_data_set.txt file in the data_sets subdirectory on the server. 

Option 2: Run Locally

You can also run Variant View locally on your own computer. 

0) Download the software package to your computer and unzip it into the desired directory. 

1) Open a terminal/command line window, and change directories to the directory where you unzipped the package.

2) Ensure that python is installed on your computer. Check by typing "python" at the command line in the terminal/console window. It is automatically installed on most versions of Mac OSX. If it is not already on your machine, you can download it from http://www.python.org/

3) Run the following command on the command line: "python -m SimpleHTTPServer 8888 &". (This command starts a python simple server to serve up the contents of the directory it is run from.)

4)  Open up a browser (Chrome, Safari, Opera, and Firefox are supported; IE is not), and point the browser to the python server with the URL "http://localhost:8888/VariantView.html"

5) You should see Variant View's display within the browser.

(III) Variant View Tool Usage

Variant View will automatically load the variant dataset you provide it, and present the contents of your file as a list of genes at the right of the display.  Selecting a gene from the list displays a visual encoding for the transcript, protein regions, and variants associated with that gene. A description of the visual encoding is provided in the legend.png file.  Alternative transcript models are available at the top of the gene display. By selecting sorting buttons above the list of genes, you can reorder the list by alphabetical order, and two other scoring metrics that rank genes in order of interesting variant distribution properties. Below the gene visual encoding view in the centre of the display there is a detailed variant view. This view shows additional variant information in tabular format.  This view is bidirectionally linked to the gene view: moving the cursor over a variant in the detailed table view will highlight the corresponding variant in the gene view and vice versa. If you already know which gene you would like to see variants in, you can type the gene name, or id, into the search box at the top of the display, and click on the search button to display it in the gene view.

For more information on Variant View, see the Variant View paper page at: http://www.cs.ubc.ca/labs/imager/tr/2013/VariantView/ ; the page contains a publication describing Variant View, and a video depicting several interesting use cases for Variant View.


(IV) Variant Data File Format Description

The provided File 6) is an example of the variant data format Variant View supports. The data file is a tab separated flat file.  Each column is a different variant attribute. Below we provide a description of what should be in each column, along with whether the column's content is absolutely required for Variant View to work (REQUIRED), or strongly suggested for meaningful analysis to take place (Suggested). The demo data file can be referenced to determine the correct format for data in each column.

Column (1) Variant ID. This field can be empty.

Column (2) Chromosome Location (Genome Coordinate)(REQUIRED)

Column (3) Ref base (Strongly Suggested)

Column (4) Var base (Strongly Suggested)

Column (5) dbSNP129 - Known harmless database 1 (Suggested)

Column (6) dbSNP135 - Known harmless database 2 (Suggested)

Column (7) dbSNP137 - Known harmless database 3 (Suggested)

Column (8) COSMIC - Known harmful/cancer database (Suggested)

Column (9) Variant Type - mutation type (Strongly Suggested)

Column (10) Amino Acid (AA) change - the protein coordinate, and the amino acid identity change at that coordinate (Strongly Suggested)

Column (11): GeneName (REQUIRED)

Column (12): Gene ID (Refseq ID Refgene ID) (REQUIRED)
