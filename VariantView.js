/*
 The first function that gets called.  We load a lot of data at the onset.
 Information from the refgene database is loaded. Information from the uniprot
 database is loaded. Variant data is also loaded at this time.
*/

// Pre-load - we do a lot of data loading on load.
window.onload = function() {
    // Load the Refgene KB:
    
    g_ref_gene = loadCSV("Knowledge_Bases/refGene.csv");
    
    // Load the Uniprot KB:
    g_uni_prot = load_uniprot_file("Knowledge_Bases/uniprotIDFTREFSEQ.dat");
    
    // This might get moved around depending on if we let users load their own files.
    // Load the variants:
    load_variant_data();
    
    // The raw variant data: variant_data_gl;
    parse_variant_data_into_gene_objects();
    
    gene_selection_display();
 
   
}


/*
 The parse_variant_data_into_gene_objects function takes a variant data set, sorts by gene id, and places all variants for a particular gene id into
 a particular gene id object. 
 */

var gene_objects;

function parse_variant_data_into_gene_objects(){
    gene_objects = [];
    
    // Do this first, and then collect the ids into gene objects
    
    var gene_id_index = ref_seq_id_index;
    var found_list = [];
    
    for(i = 0; i < variant_data_gl.length; i++){ 
        
        
        var found = found_list.filter(function(val) {
                                      return (val == variant_data_gl[i][gene_id_index]);
                                      });
        
        
        
        if(found.length == 0){
            
            var gene_id_variants = variant_data_gl.filter(function(val) {
                                                          return (val[gene_id_index] == variant_data_gl[i][gene_id_index]);
                                                          });
            
            var gene_object = new Object();
            
            
            
            gene_object.gene_name = gene_id_variants[0][gene_name_index];
            gene_object.gene_id = gene_id_variants[0][gene_id_index];
            
            
            var aa_coords = [];
            var aa_changes = [];
            var gen_coords = [];
            
            var var_objects = []; // have an object for each variant line - have to extend this to a long line.
            
            for(j = 0; j < gene_id_variants.length; j++){
                
                var current_aa_coord = gene_id_variants[j][aa_change_index]; // 70 in new file
                
                aa_changes.push(current_aa_coord);
                
                var aa_coord = Number(current_aa_coord.match(/\d+\.?\d*/g));
                
                aa_coords.push(aa_coord);
                
                var current_gen_coord = Number(gene_id_variants[j][genome_coord_index]);
                
                gen_coords.push(current_gen_coord);
                
                // The variant object stuff is to replace above:
                var var_object = new Object();
                var_object.aa_coord = aa_coord;
                var_object.aa_change = current_aa_coord;
                var_object.gen_coord = current_gen_coord;
                
                var_object.var_entry = gene_id_variants[j];
                
                var_objects.push(var_object);
                
            }
            
            // Time to do some sorting.
            var sorting_aa_indices = [];
            for(j = 0; j < aa_coords.length; j++){
                
                var list = [j,aa_coords[j]];
                
                sorting_aa_indices.push(list);
                
            }
            
            
            sorting_aa_indices.sort(function(a,b){
                                    return a[1] - b[1];
                                    });
            
            var sorted_vars = [];
            var sorted_aas = [];
            var sorted_aa_changes = [];
            var sorted_gen_coord_by_aa = [];
            
            for(j = 0; j < aa_coords.length; j++){
                var index = sorting_aa_indices[j][0];
                sorted_aas.push(aa_coords[index]);
                sorted_vars.push(gene_id_variants[index]);
                sorted_aa_changes.push(aa_changes[index]);
                sorted_gen_coord_by_aa.push(gen_coords[index]);
            }
            
            gene_id_variants = sorted_vars;
            aa_coords = sorted_aas;
            aa_changes = sorted_aa_changes;
            

            gene_object.num_variants = gene_id_variants.length;
            gene_object.variants = gene_id_variants; 
            gene_object.aa_coords = aa_coords;
            gene_object.aa_changes = aa_changes;
            gene_object.gen_coords = sorted_gen_coord_by_aa;
            
            
            //calculate the hotspot
            
            var unsorted_variants = gene_object.gen_coords;
            
            sorted_variants = unsorted_variants.sort(function(a,b){
                                                     return a - b;
                                                     });
            
            var groups = [];
            var group = 0;
            
            for(i = 0; i < sorted_variants.length; i++){
                
                if(i==0 && i==sorted_variants.length-1){
                    group = 1;
                    groups.push(group);
                }else if(i==0){
                    group = 1;
                }else if(i==sorted_variants.length-1){
                    if((sorted_variants[i] - sorted_variants[i-1]) <=20){
                        group = group+1;
                        groups.push(group);
                    }else{
                        
                        groups.push(group);
                        groups.push(1);
                    }
                }else{
                    if((sorted_variants[i] - sorted_variants[i-1]) <=20){
                        group = group+1;
                    }else{
                        groups.push(group);
                        group = 1;
                    }
                }
            }
            
            // Now get the max of these groups
            var max = 0;
            //alert("The groups: " + groups);
            for(i = 0; i < groups.length; i++){
                if(groups[i] > max){
                    max = groups[i]; 
                }
            }
            
            gene_object.hotspot_score = max;
            
            var sorted_var_object = var_objects.sort(function(a,b){
                                                     return a.gen_coord - b.gen_coord;
                                                     });
            gene_object.variant_objects = sorted_var_object;
            
            gene_objects.push(gene_object);
            
            found_list.push(gene_id_variants[0][ref_seq_id_index]);
        }
    } 
}


/*
 The showAndClearField function is called by the VariantView.html file. This
 function serves to take user gene input from the search box and query the data
 set.
 */

function showAndClearField(form){
   
    var gene_name_or_id = form.inputbox.value;
    getCurrGeneSelection(gene_name_or_id);
}
/*
 File indices:
 We set the index for the column number in the flat file where each attribute
 is located. These file column indices are for the variant data file.
*/

var patient_id_index = 0;
var genome_coord_index = 1;
var ref_base_index = 2;
var var_base_index = 3;
var dbSNP_129_index = 4;
var dbSNP_135_index = 5;
var dbSNP_137_index = 6;
var cosmic_index = 7;
var variant_type_index = 8;
var aa_change_index = 9;
var gene_name_index = 10;
var ref_seq_id_index = 11;


var filtered_variants;

function remove_current_variants(){
    
    // Remove variant lines
    remove = gene_model_display.selectAll(".variant_line_summ")
    .data([], String);
    remove.enter().append(".variant_line_summ")
    ;
    
    remove.exit().remove();
    
    // Remove variant circles
    remove = gene_model_display.selectAll(".variant_circle_summ")
    .data([], String);
    remove.enter().append(".variant_circle_summ")
    ;
    
    remove.exit().remove();
}

/* Query gene field - we use the entered gene to search (Gene or ID or both?) */

var variant_data_gl;
var variant_data_rs;
var chromosomes;
var query_gene_variants;


/*
 For queries that don't use user input from text box but gene selection buttons.
*/
 
 function getCurrGeneSelection(selection) {
    
    // Here we are expecting the "NM_..." 
    // Need to support gene name as well!
    // The refseq id is at 1; the gene name is at index 12
    g_id_or_name = 0;
    g_queryGene_name = selection;
    
    // This was the file format from gerben and linda - do we want it to be different?
    
    // First check whether it's an id or not
    if(g_queryGene_name.indexOf("NM_" || "NR_") != -1){
        // It's a gene id
        
        g_id_or_name = 0;
        
    }else{
        // otherwise it's probably not...
        
        g_id_or_name = 1;
    }
    
    // Length 65
    // Remove all NA's
    // Find the correct gene!
    
    // Need to get all gene objects - same gene names, different transcripts
    
    get_gene_object(); 

    load_ref_gene_and_uniprot(); 
     
    // Add the exome coord, exome screen coordinate, and exome number to the gene object.
    add_transcript_info_to_gene_object();
    
    data_parser(); // First time we draw variants
    redraw_model();
    gene_selection_display();
    
    // Now show the alternative transcript buttons at the top of the screen. 
     
     
    // Back to the top of the screen for a moment
    gene_model_display.append("text")
    .text("Alternative Transcripts:")
    .attr("x", 10)
    .attr("y", 20)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    // Draw them!
    var x_spacer_alt_trans = 0;
    var y_spacer_alt_trans = 0;
    for(i = 0; i < g_all_gene_objects.length; i ++){
        gene_objects[i].gene_id;
        
        var index = [i];
 
        gene_model_display.append("text")
        .text(g_all_gene_objects[i].gene_name + " (" + g_all_gene_objects[i].gene_id + ") ")
        //.text("  gene-anon" + " (trans-anon) ")
        .attr("class", "alt_trans_text" + i)
        .attr("x", 10 +200 + x_spacer_alt_trans)
        .attr("y", 20)
        .attr("font-family", "sans-serif")
        .attr("font-size", "12px")
        .attr("fill", "#666");
        
        gene_model_display
        .data(index)
        .append("rect")
        .attr("class", "alt_trans_button" + i)
        .attr("x", 10 +200 +x_spacer_alt_trans)
        .attr("y", 5 + y_spacer_alt_trans) // was 25 use this for exons
        .attr("width", 155)
        .attr("height", 20) // Was 35 use this for exons
        .attr("stroke-width", 1.5)
        .attr("stroke", "#F0F0F0")
        .attr("fill", "#F0F0F0")
        .attr("fill-opacity", 0.0)
        .attr("stroke-opacity", 1.0)
        .on("mouseover",function(d) {
            d3.select(this).style("fill", "#ddd");
            d3.select(this).style("stroke", "#000000");
            d3.select(this).style("stroke-opacity",1.0);
            })
        .on("mouseout",function(d) {
            d3.select(this).style("fill", "#F0F0F0");
            d3.select(this).style("stroke", "#F0F0F0");
            d3.select(this).style("stroke-opacity",0.50);
            })
        .on("mousedown", function() 
            {   // Draw the alternative transcripts
            d3.select(this).style("fill", "aliceblue");
            
            
            })
        .on("mouseup", function(d) 
            {  
            d3.select(this).style("fill", "#F0F0F0");
            d3.select(this).style("stroke", "#F0F0F0");
            d3.select(this).style("stroke-opacity",0.50); 
            
            getCurrGeneSelection(g_all_gene_objects[d].gene_id); 
            
            });
        
        
        x_spacer_alt_trans += 155;
        
    }
    
}

/*
 The get_gene_object function retrieves a gene object and its variants from 
 the list of gene objects and returns it.
 
 */


var g_all_gene_objects;

function get_gene_object(){
    
    // Get all of them
    g_all_gene_objects = [];
    
    if(g_id_or_name == 0){
        // then it's an id
        for(i = 0; i < gene_objects.length; i++){
            if(g_queryGene_name == gene_objects[i].gene_id){
                queried_gene_object = gene_objects[i];
            }
        }
        // Now get all that fall under a gene name:
        var gene_name = queried_gene_object.gene_name;
        
        for(i = 0; i < gene_objects.length; i++){
            if(gene_name == gene_objects[i].gene_name){
                g_all_gene_objects.push(gene_objects[i]);
            }
        }
        
        
    }else{
        // else it's a name
        for(i = 0; i < gene_objects.length; i++){
            if(g_queryGene_name == gene_objects[i].gene_name){
                
                g_all_gene_objects.push(gene_objects[i]);
                queried_gene_object = gene_objects[i];
            } 
        }
    }
}

/*
 The add_transcript_info_to_gene_object retrieves information from
 transcript model for display.
 
 */

function add_transcript_info_to_gene_object(){
    var transcript_model = new Object();
    
    refgene_entry_parsing(); // Takes care of where genome coordinates land in exons and parsing refgene entry.
    get_exome_position();
    
    queried_gene_object.sense = transcript_sense;
    
    // Correct the exome coordinates for the screen display (left to right with the correct sense)
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i ++){
        var exon_coord = 0;
        
        if(queried_gene_object.sense == "-"){
            exon_coord = exon_model_length - queried_gene_object.variant_objects[i].exome_coord;
        }else{
            exon_coord = queried_gene_object.variant_objects[i].exome_coord;  
        }
        
        queried_gene_object.variant_objects[i].exome_coord = exon_coord;
    }
}


/*
This load_variant_data function loads the tab-separated text file that contains
your variant data. To change the name of the file, edit the string value of the
global variable "variant_file_name"

*/

var variant_file_name = "demo_data_set.txt"
var whitespace_removed_gl_var;
function load_variant_data(){
    var gl_file = variant_file_name;
    
    variant_data_gl = load_gl_("data_files/" + gl_file);
    
    // remove the first line it's a header
    variant_data_gl.splice(0,1);
    
    return 1;
}

/*
 Load function associated with the above load_variant_data function.
 
 */

function load_gl_(file) {
    if (window.XMLHttpRequest) {
        // IE7+, Firefox, Chrome, Opera, Safari
        var request = new XMLHttpRequest();
    }else{
        // code for IE6, IE5
        var request = new ActiveXObject('Microsoft.XMLHTTP');
    }
    // load
    request.open('GET', file, false);
    request.send();
    
    return parse_gl_(request.responseText);
}

/*
 The parse_gl_ function splits the data into an array of lines. The flat variant
 text file becomes a 2 dimensional array.
 
 */

function parse_gl_(data){ 
    //replace UNIX new lines
    data = data.replace (/\r\n/g, "\n");
    //replace MAC new lines
    data = data.replace (/\r/g, "\n");
    //split into rows
    var rows = data.split("\n");
    // create array which will hold our data:
    var gl_data = [];
    
    // loop through all rows
    for (var i = 0; i < rows.length; i++){
        // this line helps to skip empty rows
        if (rows[i]) {                    
            // our columns are separated by comma
            // Each line is: "Sample ID, Pos ...."
            var column = rows[i].split("\t");  
            
            
            
            gl_data.push(column);
        }
    }
    
    return gl_data;
}

// Loads both refgene transcript model data and protein model data
function load_ref_gene_and_uniprot(){
    var success_ref_gene = load_ref_gene();
    if(success_ref_gene){
        // Success! We're going to query for the protein model now!
        
        var success = load_uni_prot();
        if(success == 1){
            //alert("Successfully loaded transcript and protein model");
            return 1;
        }else{
            alert("No UniProt record available");
            return 0;
        }
    }else{
        return 0;
        // Do nothing
    }
    
    return 0;
}

var g_ref_gene;
var g_ref_gene_entry;

var g_id_or_name;

/*
 The load_ref_gene function uses a gene id to retrieve information from refseq 
 entry; these can include coordinates for exons.
 
 */

function load_ref_gene(){
    
    g_ref_gene_entry = retrieve_refGene_entry_ID(queried_gene_object.gene_id);
    
    if(g_ref_gene_entry == 0){
        alert("Could Not Find RefGene ID: " + g_queryGene_name + " - Please Try Again");
        return 0;
    }else{
       
        
        return 1;
    }

    return 1;
}

function retrieve_refGene_entry_ID(p_id){
    var matching_refgene_transcripts = [];
    if(g_id_or_name == 0){
        
        for(i = 0; i < g_ref_gene.length; i++){
            if((p_id == g_ref_gene[i][1]) && (g_ref_gene[i][2].indexOf('_') == -1)){
                //alert(g_ref_gene[i][2]);
                matching_refgene_transcripts.push(g_ref_gene[i]);
                //return g_ref_gene[i];
            }else{
                // Have not found it yet
            }
        }
        
    }else{
        
        for(i = 0; i < g_ref_gene.length; i++){
            // HACK!
            if(g_ref_gene[i].indexOf(p_id) != -1){
                matching_refgene_transcripts.push(g_ref_gene[i]);
                //return g_ref_gene[i];
            }else{
                // Have not found it yet
            }
        } 
    }
    
    for(j = 0; j < g_all_gene_objects.length; j++){
        
        
        for(i = 0; i < g_ref_gene.length; i++){
            if((g_all_gene_objects[j].gene_id == g_ref_gene[i][1]) && (g_ref_gene[i][2].indexOf('_') == -1)){
                
                g_all_gene_objects[j].transcript_info = g_ref_gene[i];
                // alert(g_all_gene_objects[j].transcript_info);
                matching_refgene_transcripts.push(g_ref_gene[i]);
                
            }else{
                // Have not found it yet
            }
        }
        
    }
    
    //alert(g_all_gene_objects.length);
    
    if(matching_refgene_transcripts.length == 0){
        return 0;
    }else{
        return matching_refgene_transcripts[0];
    }
    
}

var g_uni_prot;
var g_uni_prot_entry;

/*
 This function retrieves associated uniprot information for the refseq id and 
 transcript.
 
 */

function load_uni_prot(){
    //alert("Loading UniProt Knowledge Base - Click 'Okay'");
    //g_uni_prot = load_uniprot_file("Knowledge_Bases/uniprotIDFTREFSEQ.dat");
    g_uni_prot_entry = retrieve_uniprot_entry(queried_gene_object.gene_id);
    if(g_uni_prot_entry == 0){
        alert("No UniProt Record for: " + g_queryGene_name);
        return 0;
    }else{
        //alert("Found Protein Model");
        // Wipe the rest of uniprot from memory we don't need it
        // g_uni_prot = [];
        
        return 1;
    }
    return 1;
}

/* 
 The retrieve_uniprot_entry function is a helper function for the above uniprot
 retrieval function.
 */

function retrieve_uniprot_entry(p_refSeqId){
    var uni_prot_length = g_uni_prot.length;
    
    // Need a "."!
    
    p_refSeqId += ".";
    
    for(i = 0; i < uni_prot_length; i++){
        if(g_uni_prot[i].indexOf(p_refSeqId) != -1){
            // We need to keep moving up to the next ID - this will ensure we
            // get the whole record
            g_uni_prot_entry = []
            var index = 0;
            // Remember to put "E" at end of file!
            while(g_uni_prot[i + index][0] != "I"){
                g_uni_prot_entry.push(g_uni_prot[i + index]);
                index = index + 1;
            }
            return g_uni_prot_entry;
        }else{
            // Have not found it yet
        }
    }
    return 0;
}


function load_uniprot_file(file){
    if (window.XMLHttpRequest) {
        // IE7+, Firefox, Chrome, Opera, Safari
        var request = new XMLHttpRequest();
    }
    else {
        // code for IE6, IE5
        var request = new ActiveXObject('Microsoft.XMLHTTP');
    }
    // load
    request.open('GET', file, false);
    request.send();
    
    return parse_uniprot(request.responseText);
}

function parse_uniprot(data){ 
    //replace UNIX new lines
    data = data.replace (/\r\n/g, "\n");
    //replace MAC new lines
    data = data.replace (/\r/g, "\n");
    //split into rows
    var rows = data.split("\n");
    // create array which will hold our data:
    dataProvider = [];
    
    // loop through all rows
    for (var i = 0; i < rows.length; i++){
        // this line helps to skip empty rows
        if (rows[i]) {                    
            
            dataProvider.push(rows[i]);
        }
    }
    
    return dataProvider;
}

var g_filtered_variant_list;




var gene_display_x_offset = 110;
var gene_display_y_offset = 0;
var amino_acid_detail_x_offset = 10;

/*
Screen dimensions. These could be changed based on user preferences.
 
 */
var screen_width_protein_model = 675;

var screen_width_gene_model = 675;

var recurrence_hist_width = 675;


var middle_pane_top_y_spacer = 50;

// Calls the drawing methods
function redraw_model(){
    
    
    
    remove_all();
    // Background outlines, text etc...
    drawPermanentStructures();
    
    drawExonSummaryView();
    
    drawProteinSummaryView(); 
    
    generate_variant_list();
    
}

// Create a blank canvas to draw on again

function remove_all(){
    // removal
    
    var remove = gene_model_display.selectAll("rect")
    .data([], String);
    remove.enter().append("rect")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll("text")
    .data([], String);
    remove.enter().append("text")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll("line")
    .data([], String);
    remove.enter().append("line")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll("circle")
    .data([], String);
    remove.enter().append("circle")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll("path")
    .data([], String);
    remove.enter().append("path")
    ;
    
    remove.exit().remove();
    
}

// global to hold the name of the entered data set
var dataset_name;
var currentDataset;
var g_queryGene_name;
var gene_model;

// Holds Values of the Sliders
var quality_filter;
var depth_filter;

var rawData, // the raw data from the csv file
drawingData, // Data we don't want to display (dirty or it's type is unchecked)
currentDataset, // name of the current data set. Used to track when the data set changes.
gene_model_display,
compound_list_display,
detailed_list_display,
gene_model_display; // the visualization selection



// Holding our data!
var variantData;
var refGeneData;
var uniprotData;


function loadCSV(file) {
    if (window.XMLHttpRequest) {
        // IE7+, Firefox, Chrome, Opera, Safari
        var request = new XMLHttpRequest();
    }
    else {
        // code for IE6, IE5
        var request = new ActiveXObject('Microsoft.XMLHTTP');
    }
    // load
    request.open('GET', file, false);
    request.send();
    
    return parseCSV(request.responseText);
}

function parseCSV(data){ 
    //replace UNIX new lines
    data = data.replace (/\r\n/g, "\n");
    //replace MAC new lines
    data = data.replace (/\r/g, "\n");
    //split into rows
    var rows = data.split("\n");
    // create array which will hold our data:
    dataProvider = [];
    
    // loop through all rows
    for (var i = 0; i < rows.length; i++){
        // this line helps to skip empty rows
        if (rows[i]) {                    
            // our columns are separated by comma
            // Each line is: "Sample ID, Pos ...."
            var column = rows[i].split(",");  
            
            dataProvider.push(column);
        }
    }
    
    return dataProvider;
}



var transcriptStart;
var transcriptEnd;
var codingRegStart;
var codingRegEnd;
var numberOfExons;
var length_of_splice_site;
var total_number_splices;
var exonStarts;
var exonEnds;
var geneName;

// A lot of alt transcript stuff
var numberAltTranscripts;

// There should probably be some objec that holds all the meta data for the
// the alt trans and protein models.

var exon_splicesite_screencoords,exonLengthsScreen,bp_length_model;

// Function to filter and update the visualization based on controls
function filterDrawUpdate(){
    
}

var signals;
var whole_chains;
var domains;
var regions;
var compbiases;
var act_sites;
var metals;
var bindings;
var lipids;
var carbohyds;
var topo_domains;
var transmembranes;
var disulfides;
var zinc_fingers;
var np_bind;
var mod_res;

// Annotation for protein

var domain_info;
var topo_info;
var trans_info;
var aa_chain_info;
var region_info;
var comp_info;
var act_info;
var metal_info;
var binding_info;
var lipid_info;
var carb_info;
var disulf_info;
var signal_info;
var zingfing_info;
var npbind_info;
var modres_info;

var length_gene_model_splices_included;

var prot_annot_info;

// Called when we load a new data set - typically all new transcripts
// populate the model
var exon_lengths;

function data_parser(){
    // Reset all the globals
    
    signals = [];
    whole_chains = [];
    domains = [];
    regions = [];
    compbiases = [];
    act_sites = [];
    metals = [];
    bindings = [];
    lipids = [];
    carbohyds = [];
    topo_domains = [];
    transmembranes = [];
    disulfides = [];
    zinc_fingers = [];
    np_bind = [];
    mod_res = [];
    
    // Some annotations
    domain_info = [];
    topo_info = [];
    trans_info = [];
    aa_chain_info = [];
    region_info = [];
    comp_info = [];
    act_info = [];
    metal_info = [];
    binding_info = [];
    lipid_info = [];
    carb_info = [];
    disulf_info = [];
    signal_info = [];
    zingfing_info = [];
    npbind_info = [];
    modres_info = [];
    
    
    for(i = 0; i < g_uni_prot_entry.length; i++){
        if(g_uni_prot_entry[i].indexOf("CHAIN") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            
            aa_chain_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            whole_chains.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("DOMAIN") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            domain_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            domains.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("REGION") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            region_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            regions.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("COMPBIAS") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            comp_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            compbiases.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("ACT_SITE") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            act_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            act_sites.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("METAL") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            metal_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            metals.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("BINDING") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            binding_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            bindings.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("LIPID") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            lipid_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            lipids.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("CARBOHYD") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            carb_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            carbohyds.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("TOPO_DOM") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            topo_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            topo_domains.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("TRANSMEM") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            
            trans_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            
            transmembranes.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("DISULFID") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            var info = g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length);
            if(info == ""){
                info = "Disulfide Bond";
            }
            disulf_info.push(info);
            
            disulfides.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("SIGNAL") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            signal_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            
            signals.push(interval);
        }else if (g_uni_prot_entry[i].indexOf("ZN_FING") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            zingfing_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            
            zinc_fingers.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("NP_BIND") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            npbind_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            
            np_bind.push(interval);
        }else if(g_uni_prot_entry[i].indexOf("MOD_RES") != -1){
            var interval = get_uniprot_FT_interval(g_uni_prot_entry[i]);
            modres_info.push(g_uni_prot_entry[i].substring(30,g_uni_prot_entry[i].length));
            
            mod_res.push(interval);
        }else{
            
        }
        
    }
    protein_annotations = [];
    
    // Now put these lists of protein annotations in a large list for later
    protein_annotations.push(whole_chains);
    protein_annotations.push(signals);
    protein_annotations.push(domains);
    protein_annotations.push(regions);
    protein_annotations.push(compbiases);
    protein_annotations.push(topo_domains);
    protein_annotations.push(transmembranes);
    protein_annotations.push(zinc_fingers);
    protein_annotations.push(act_sites);
    protein_annotations.push(np_bind);
    protein_annotations.push(metals);
    protein_annotations.push(bindings);
    protein_annotations.push(mod_res);
    protein_annotations.push(lipids);
    protein_annotations.push(carbohyds);
    protein_annotations.push(disulfides);
    
    prot_annot_info = [];
    
    prot_annot_info.push(aa_chain_info);
    prot_annot_info.push(signal_info);
    
    prot_annot_info.push(domain_info);
    prot_annot_info.push(region_info);
    prot_annot_info.push(comp_info);
    prot_annot_info.push(topo_info);
    prot_annot_info.push(trans_info);
    
    prot_annot_info.push(zingfing_info);
    prot_annot_info.push(act_info);
    prot_annot_info.push(npbind_info);
    prot_annot_info.push(metal_info);
    prot_annot_info.push(binding_info);
    prot_annot_info.push(modres_info);
    
    prot_annot_info.push(lipid_info);
    prot_annot_info.push(carb_info);
    prot_annot_info.push(disulf_info);
    
}

var transcript_sense;
function refgene_entry_parsing(){
    transcript_sense = g_ref_gene_entry[3];
    transcriptStart = Number(g_ref_gene_entry[4]);
    transcriptEnd = Number(g_ref_gene_entry[5]);
    codingRegStart = Number(g_ref_gene_entry[6]);
    codingRegEnd = Number(g_ref_gene_entry[7]);
    numberOfExons = Number(g_ref_gene_entry[8]);
        
    //Exon Starts
    var index = 1 + 9;
    exonStarts = [];
    
    exonStarts.push(Number(g_ref_gene_entry[9].replace('"',''))); // Get first one
    
    while(g_ref_gene_entry[index].indexOf('"') == -1){
        exonStarts.push(Number(g_ref_gene_entry[index]));    
        index = index + 1;
    }
    
    exonStarts.push(Number(g_ref_gene_entry[index].replace('"',''))); // Get last one
    exonStarts.pop();
    
    // Now for the ends
    index = index + 1;
    exonEnds = [];
    exonEnds.push(Number(g_ref_gene_entry[index].replace('"',''))); // Get first one
    index = index + 1;
    while(g_ref_gene_entry[index].indexOf('"') == -1){
        exonEnds.push(Number(g_ref_gene_entry[index]));    
        index = index + 1;
    }
    
    exonEnds.push(Number(g_ref_gene_entry[index].replace('"',''))); // Get last
    
    exonEnds.pop();
    
    // Now to get anything else, the index will start at...
    index = index + 2;
    
    geneName = g_ref_gene_entry[index];
    
    exon_lengths = []
    exon_model_length = 0;
    
    var starts = exonStarts;
    
    var ends = exonEnds;
    
    for(i = 0; i < numberOfExons; i++){
        if(exonStarts[i] < codingRegStart){
            if(exonEnds[i] < codingRegStart){
                // Remove this exon
                starts.splice(i,1);
                ends.splice(i,1);
            }else{
                // Clip this exon.
                starts[i] = codingRegStart;
            }
        }
        
        if(exonEnds[i] > codingRegEnd){
            if(exonStarts[i] > codingRegEnd){
                // Remove this exon
                starts.splice(i,1);
                ends.splice(i,1);
            }else{
                // Clip this exon.
                ends[i] = codingRegEnd;
            }
        }
        
    }
    
    exonStarts = starts;
    exonEnds = ends;
    
    numberOfExons = exonEnds.length
    
    for(i = 0; i < numberOfExons; i++){
        exon_model_length += exonEnds[i] - exonStarts[i];
        exon_lengths.push(exonEnds[i] - exonStarts[i]);
    }
    
    // Clip exons that are not in the coding region.
    
    // UPDATE do this for multiple transcripts
    
    for(j = 0; j < g_all_gene_objects.length; j++){
        var transcript_entry = g_all_gene_objects[j].transcript_info;
        
        g_all_gene_objects[j].transcript_sense = transcript_entry[3];
        g_all_gene_objects[j].transcriptStart = Number(transcript_entry[4]);
        g_all_gene_objects[j].transcriptEnd = Number(transcript_entry[5]);
        g_all_gene_objects[j].codingRegStart = Number(transcript_entry[6]);
        g_all_gene_objects[j].codingRegEnd = Number(transcript_entry[7]);
        
        
        
        var numberOfExons_ = Number(transcript_entry[8]);
                
        //Exon Starts
        var index = 1 + 9;
        var exonStarts_ = [];
        
        exonStarts_.push(Number(transcript_entry[9].replace('"',''))); // Get first one
        
        while(transcript_entry[index].indexOf('"') == -1){
            exonStarts_.push(Number(transcript_entry[index]));    
            index = index + 1;
        }
        
        exonStarts_.push(Number(transcript_entry[index].replace('"',''))); // Get last one
        exonStarts_.pop();
        
        // Now for the ends
        index = index + 1;
        var exonEnds_ = [];
        exonEnds_.push(Number(transcript_entry[index].replace('"',''))); // Get first one
        index = index + 1;
        while(transcript_entry[index].indexOf('"') == -1){
            exonEnds_.push(Number(transcript_entry[index]));    
            index = index + 1;
        }
        
        exonEnds_.push(Number(transcript_entry[index].replace('"',''))); // Get last
        
        exonEnds_.pop();
        
        // Now to get anything else, the index will start at...
        index = index + 2;
        
        //geneName = transcript_entry[index];
        
        var exon_lengths_ = []
        var exon_model_length_ = 0;
        
        var starts = exonStarts_;
        
        var ends = exonEnds_;
        
        for(i = 0; i < numberOfExons_; i++){
            if(exonStarts_[i] < codingRegStart){
                if(exonEnds_[i] < codingRegStart){
                    // Remove this exon
                    starts.splice(i,1);
                    ends.splice(i,1);
                }else{
                    // Clip this exon.
                    starts[i] = codingRegStart;
                }
            }
            
            if(exonEnds_[i] > codingRegEnd){
                if(exonStarts_[i] > codingRegEnd){
                    // Remove this exon
                    starts.splice(i,1);
                    ends.splice(i,1);
                }else{
                    // Clip this exon.
                    ends[i] = codingRegEnd;
                }
            }
            
        }
        
        exonStarts_ = starts;
        exonEnds_ = ends;
        
        numberOfExons_ = exonEnds_.length
        
        for(i = 0; i < numberOfExons_; i++){
            exon_model_length_ += exonEnds_[i] - exonStarts_[i];
            exon_lengths_.push(exonEnds_[i] - exonStarts_[i]);
        }
        
        g_all_gene_objects[j].exonStarts = exonStarts_;
        g_all_gene_objects[j].exonEnds = exonEnds_;
        g_all_gene_objects[j].exonLengths = exon_lengths_;
        g_all_gene_objects[j].numberOfExons = numberOfExons_;
        
    }
    
}

/*
 The get_exome_position converts a genome coordinate to an exome coordinate. This conversion requires knowledge of the genome coordinate AND the associated transcript coordinates for the gene.
 
 */

function get_exome_position(){
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i++){
        
        var gaps = 0;
        for(j = 0; j < numberOfExons; j++){
            
            if(j == 0){ // First exon
                
                if((exonStarts[j] <= queried_gene_object.variant_objects[i].gen_coord) && (exonEnds[j] >= queried_gene_object.variant_objects[i].gen_coord)){
                    //alert("One");
                    var coord = queried_gene_object.variant_objects[i].gen_coord - exonStarts[0];
                    //alert(coord);
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }else if(exonStarts[j] > queried_gene_object.variant_objects[i].gen_coord){
                    //alert("Two");
                    var coord = 0;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }
            }else if(j > 0 && j < numberOfExons-1){ // Middle Exons
                gaps += exonStarts[j] - exonEnds[j-1];
                if(exonStarts[j] <= queried_gene_object.variant_objects[i].gen_coord && exonEnds[j] >= queried_gene_object.variant_objects[i].gen_coord){
                    
                    
                    var coord = queried_gene_object.variant_objects[i].gen_coord - exonStarts[0] - gaps;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }else if((exonEnds[j-1] < queried_gene_object.variant_objects[i].gen_coord) && (exonStarts[j] > queried_gene_object.variant_objects[i].gen_coord)){
                    // alert("Four");
                    
                    var coord = exonStarts[j] - exonStarts[0] - gaps;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }
            }else if(j == (numberOfExons - 1)){ // Last exon
                gaps += exonStarts[j] - exonEnds[j-1];
                if((exonStarts[j] <= queried_gene_object.variant_objects[i].gen_coord) && (exonEnds[j] >= queried_gene_object.variant_objects[i].gen_coord)){
                    //alert("Five");
                    
                    var coord = queried_gene_object.variant_objects[i].gen_coord - exonStarts[0] - gaps;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }else if((exonEnds[j-1] < queried_gene_object.variant_objects[i].gen_coord) && (exonStarts[j] > queried_gene_object.variant_objects[i].gen_coord)){
                    
                    //alert("Six");
                    
                    var coord = exonStarts[j] - exonStarts[0] - gaps;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }else if(queried_gene_object.variant_objects[i].gen_coord > exonEnds[numberOfExons -1]){
                    //alert("Seven");
                    
                    var coord = exon_model_length;
                    
                    queried_gene_object.variant_objects[i].exome_coord = coord;
                }
            }
        }
    }
}

var queried_gene_object;

var variants_with_exon_coords;

var filtered_samples_by_exon_total_list;
// Grab the intervals from a uniprot line
function get_uniprot_FT_interval(p_uniprot_line){
    var split_string = p_uniprot_line.split(" ");
    var removed_white = split_string.filter(function(val) {
                                            return val != "";
                                            });
    var interval = [];
    //Force these to numbers
    interval.push(Number(removed_white[2]));
    interval.push(Number(removed_white[3]));
    return interval;
}

function get_protein_annotation_info(p_uniprot_line){
    var split_string = p_uniprot_line.split(" ");
    var removed_white = split_string.filter(function(val) {
                                            return val != "";
                                            });
    
    
    // alert(removed_white);   
}

var number_patients;
var patient_ids;
function get_number_patients(){
    var list = [];
    
    for(i = 0; i < variants_with_exon_coords.length; i++){
        list.push(variants_with_exon_coords[i][patient_id_index]);
    }
    var uniqueNames = [];
    $.each(list, function(i, el){
           if($.inArray(el, uniqueNames) === -1) uniqueNames.push(el);
           });
    
    list = [];
    number_patients = uniqueNames.length;
    patient_ids = uniqueNames;
    uniqueArray = [];
    display_number_patients();
}

function display_number_patients(){
    controls_panel.append("text")
    
    .text(number_patients)
    .attr("class","sample_text")
    // .attr("text-anchor", "middle")
    .attr("x", 70)
    .attr("y", 168);
}

var controls_panel;

// The below function can go! TOGO

function init_gene_model () {
    // This is ensures you can't keep clicking and overlaying
    // multiple copies of the alternative transcript view.
    alts_shown = 0; // Important switch, don't delete
    
    gene_model_display = d3.select("#compound-visualization");
    
    sortable_list_display = d3.select("#sorting-options");
    
    gene_selection_list = d3.select("#gene-selection-list");
    
    variant_table = d3.select("#variant-table");
    
    controls_panel = d3.select("#quality-controls");
    
    
    controls_panel.append("text")
    .text("Exome Forward Reads")
    // .attr("text-anchor", "middle")
    .attr("x", 36)
    .attr("y", 13)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Exome Reverse Reads")
    // .attr("text-anchor", "middle")
    .attr("x", 286)
    .attr("y", 13)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("(Forward + Reverse)/(Total)")
    // .attr("text-anchor", "middle")
    .attr("x", 534)
    .attr("y", 13)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Non-Synonymous")
    // .attr("text-anchor", "middle")
    .attr("x", 30)
    .attr("y", 100)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Frame Shift")
    // .attr("text-anchor", "middle")
    .attr("x", 190)
    .attr("y", 100)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Codon Deletion")
    // .attr("text-anchor", "middle")
    .attr("x", 305)
    .attr("y", 100)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Codon Insertion")
    // .attr("text-anchor", "middle")
    .attr("x", 440)
    .attr("y", 100)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Truncating/Nonsense")
    // .attr("text-anchor", "middle")
    .attr("x", 565)
    .attr("y", 100)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
   
    controls_panel.append("text")
    .text("Variant Recurrence at Transcript Positions")
    // .attr("text-anchor", "middle")
    .attr("x", gene_display_x_offset)
    .attr("y", 148)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Samples:")
    // .attr("text-anchor", "middle")
    .attr("x", 2)
    .attr("y", 168)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Variants:")
    // .attr("text-anchor", "middle")
    .attr("x", 2)
    .attr("y", 198)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("Unique:")
    // .attr("text-anchor", "middle")
    .attr("x", 2)
    .attr("y", 228)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    
    controls_panel.append("text")
    .text("100%")
    // .attr("text-anchor", "middle")
    .attr("x", 218 + screen_width_gene_model)
    .attr("y", 150)
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .text("0%")
    // .attr("text-anchor", "middle")
    .attr("x", 228 + screen_width_gene_model)
    .attr("y", 240)
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "#666");
    
    controls_panel.append("rect")
    .attr("x", 36) // 1
    .attr("y", 16) // 1
    .attr("width",200 + 9.1111)
    .attr("height",40)
    .attr("stroke", "#666")
    .attr("fill","#F0F0F0")
    .attr("stroke-width", "1");
    
    controls_panel.append("rect")
    .attr("x", 286)
    .attr("y", 16)
    .attr("width",200 + 9.1111)
    .attr("height",40)
    .attr("stroke", "#666")
    .attr("fill","#F0F0F0")
    .attr("stroke-width", "1");
    
    controls_panel.append("rect")
    .attr("x", 536)
    .attr("y", 16)
    .attr("width",200 + 9.1111)
    .attr("height",40)
    .attr("stroke", "#666")
    .attr("fill","#F0F0F0")
    .attr("stroke-width", "1");
    
    
    // Recurrence Box
    controls_panel.append("rect")
    .attr("x", gene_display_x_offset)
    .attr("y", 152)
    .attr("width",screen_width_gene_model)
    .attr("height",90)
    .attr("stroke", "#666")
    .attr("fill","#F0F0F0")
    .attr("stroke-width", "1");
    
    // Non Syn Check Box
    
    nonsynonymous_filter_val = 1;
    
    controls_panel.append("rect")
    .attr("x", 3)
    .attr("y", 87)
    .attr("width",15)
    .attr("height",15)
    .attr("stroke", "#666")
    .attr("fill","#666")
    .attr("stroke-width", "1")
    .on("mousedown", function() 
        {   
        if(nonsynonymous_filter_val == 1){
        d3.select(this).style("fill", "#F0F0F0");
        nonsynonymous_filter_val = 0;
        nonsynonymous_filter(); // effect the change
        }else{
        d3.select(this).style("fill", "#666");
        nonsynonymous_filter_val = 1;
        nonsynonymous_filter(); // effect the change
        }
        
        });
    
    // Frame Shift Check Box
    
    frame_shift_val = 1;
    
    controls_panel.append("rect")
    .attr("x", 166)
    .attr("y", 87)
    .attr("width",15)
    .attr("height",15)
    .attr("stroke", "#666")
    .attr("fill","#666")
    .attr("stroke-width", "1")
    .on("mousedown", function() 
        {   
        if(frame_shift_val == 1){
        d3.select(this).style("fill", "#F0F0F0");
        frame_shift_val = 0;
        
        frame_shift_filter(); // effect the change
        }else{
        d3.select(this).style("fill", "#666");
        frame_shift_val = 1;
        
        frame_shift_filter(); // effect the change
        }
        
        });
    
    // Codon Deletion Check box
    
    codon_deletion_filter_val = 1;
    
    controls_panel.append("rect")
    .attr("x", 283)
    .attr("y", 87)
    .attr("width",15)
    .attr("height",15)
    .attr("stroke", "#666")
    .attr("fill","#666")
    .attr("stroke-width", "1")
    .on("mousedown", function() 
        {   // Draw the alternative transcripts
        if(codon_deletion_filter_val == 1){
        d3.select(this).style("fill", "#F0F0F0");
        codon_deletion_filter_val = 0;
        codon_deletion_filter(); // effect the change
        }else{
        d3.select(this).style("fill", "#666");
        codon_deletion_filter_val = 1;
        codon_deletion_filter(); // effect the change
        }
        
        });
    
    
    // Codon Insertion Check box
    codon_insertion_filter_val = 1;
    
    controls_panel.append("rect")
    .attr("x", 415)
    .attr("y", 87)
    .attr("width",15)
    .attr("height",15)
    .attr("stroke", "#666")
    .attr("fill","#666")
    .attr("stroke-width", "1")
    .on("mousedown", function() 
        {   // Draw the alternative transcripts
        if(codon_insertion_filter_val == 1){
        d3.select(this).style("fill", "#F0F0F0");
        codon_insertion_filter_val = 0;
        // codon_deletion_filter(); // effect the change
        }else{
        d3.select(this).style("fill", "#666");
        codon_insertion_filter_val = 1;
        // codon_deletion_filter(); // effect the change
        }
        
        });
    
    // Truncating/Nonsense
    
    truncating_filter_val = 1;
    
    controls_panel.append("rect")
    .attr("x", 548)
    .attr("y", 87)
    .attr("width",15)
    .attr("height",15)
    .attr("stroke", "#666")
    .attr("fill","#666")
    .attr("stroke-width", "1")
    .on("mousedown", function() 
        {   // Draw the alternative transcripts
        if(truncating_filter_val == 1){
        d3.select(this).style("fill", "#F0F0F0");
        truncating_filter_val = 0;
        truncating_filter(); // effect the change
        }else{
        d3.select(this).style("fill", "#666");
        truncating_filter_val = 1;
        truncating_filter(); // effect the change
        }
        
        });
    
    
} 
var truncating_filter_val;
var codon_insertion_filter_val;
var frame_shift_val;

// Need a method to remove all this stuff!

var forward_HQ_max;
var reverse_HQ_max;
var ratio_HQ_max;


var g_max_height_bar_genome_view; // needs to be a global to access within callbck

var exon_model_length;

var s_Forward_HQ; // global brushing value on the forward widget
var s_Reverse_HQ; // global brushing value on the reverse widget
var s_Ratio_HQ; // global brushing value on the ration widget


// This can go too! TOGO
function init_control_panel(){
    
    var width = 300,
    height = 300;
    
    // Defualt threshold settings:
    
    // Pulled from forward bar
    
    g_max_height_bar_genome_view = 0;
    
    var num_bins = 200;
    
    var interval = Math.ceil(exon_model_length/num_bins);
    var interval_start = 0;
    var interval_end = 0+interval;
    var binned_counts_hist = [];
    var total_vars_possible = variants_with_exon_coords.length;
    for(i = 0; i < num_bins; i++){
        // Filtering on position
        var filts = variants_with_exon_coords.filter(function(val) {
                                                     return ((Number(val[64]) <=interval_end)&&(Number(val[64]) >= interval_start));
                                                     });
        
        
        binned_counts_hist.push(filts.length);
        interval_start = interval_end;
        interval_end = interval_end + interval;
        
    }
    
    g_max_height_bar_genome_view = Math.max.apply( Math, binned_counts_hist );
    
    // pulled from forward bar.
    
    var max = Number(variants_with_exon_coords[0][34]);
    
    for(i = 0; i < variants_with_exon_coords.length; i++){
        if((Number(variants_with_exon_coords[i][34]) > max)){
            max = Number(variants_with_exon_coords[i][34]);
        }else{
            
        }
    }
    
    forward_HQ_max = max;
    
    // Set default
    // They use 3 reads
    var forward_threshold = (3/max)
    
    // Set default threshold values - need to do this appropriately!
    s_Forward_HQ = [];
    
    s_Forward_HQ.push(forward_threshold);
    s_Forward_HQ.push(1.0);
    
    
    forward_control_bar();
    
    var max = Number(variants_with_exon_coords[0][35]);
    
    for(i = 0; i < variants_with_exon_coords.length; i++){
        if((Number(variants_with_exon_coords[i][35]) > max)){
            max = Number(variants_with_exon_coords[i][35]);
        }else{
            
        }
    }
    
    reverse_HQ_max = max;
    
    var reverse_threshold = (3/max)
    
    s_Reverse_HQ = [];
    s_Reverse_HQ.push(reverse_threshold);
    s_Reverse_HQ.push(1.0);
    
    reverse_control_bar();
    
    for(i = 0; i < variants_with_exon_coords.length; i++){
        var value1 = Number(variants_with_exon_coords[i][32]);
        var value2 = Number(variants_with_exon_coords[i][33]);
        var value3 = Number(variants_with_exon_coords[i][34]);
        var value4 = Number(variants_with_exon_coords[i][35]);
        
        
        if(isNaN(value1)|| isNaN(value2) || isNaN(value3) || isNaN(value4)){
            variants_with_exon_coords[i].push(0);
        }else{
            
            var ratio_value = (value3+value4)/(value3+value4 + value1+value2);
            variants_with_exon_coords[i].push(ratio_value*100);
            
        }
    }
    
    
    // 65 is now the value of the ratio
    var max = Number(variants_with_exon_coords[0][65]);
    
    for(i = 0; i < variants_with_exon_coords.length; i++){
        if((Number(variants_with_exon_coords[i][65]) > max)){
            max = Number(variants_with_exon_coords[i][65]);
        }else{
            
        }
    }
    
    ratio_HQ_max = max;
    
    var ratio_threshold = (10/max)
    
    s_Ratio_HQ = [];
    s_Ratio_HQ.push(ratio_threshold);
    s_Ratio_HQ.push(1.0);
    
    ratio_control_bar();
    
    recurrence_histogram_bar(true);     // Initial call
    
    remove_matrix_canvas();
    setup_variant_effect_matrix();
    
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    setup_variant_effect_matrix();
    remove_matrix_marks();
    draw_matrix_marks();
}



// Last call before we wait for interaction related call-backs

function setup_variant_effect_matrix(){
    
    gene_model_display.append("rect")
    .attr("class","matrix_canvas")
    .attr("x", gene_display_x_offset)
    .attr("y", 370)
    .attr("width",screen_width_gene_model)
    .attr("height", number_patients*2)
    .attr("stroke", "#666")
    .attr("fill","#F0F0F0")
    .attr("stroke-width", "1");
}

// Update the effect matrix with marks
function draw_matrix_marks(){
    // for each sample create a mark for the type of mutation
    
    for(i = 0; i < patient_ids.length; i++){
        
        
        var patient_data = filtered_variants.filter(function(val) {
                                                    return val[0] == patient_ids[i];
                                                    });
        
        
        for(j = 0; j < patient_data.length; j++){
            
            var variant_pos_screen = (Number(patient_data[j][64])/exon_model_length)*screen_width_gene_model;
            
            // Set some colour!
            var colour = "#666";
            //var effects = patient_data[j][11].split(",");
            //alert(effects);
            
            if(patient_data[j][11].indexOf("NON_SYNONYMOUS_CODING(") != -1){
                colour = "#FF0000";
            }else{
                
            }
            
            
            gene_model_display.append("rect")
            .attr("class","matrix_mark")
            .attr("x", gene_display_x_offset + variant_pos_screen)
            .attr("y", 370 + i*2)
            .attr("width",1)
            .attr("height", 1)
            .attr("stroke", colour)
            .attr("fill",colour)
            
            
        }
    }
}

function remove_matrix_canvas(){
    remove = gene_model_display.selectAll(".matrix_canvas")
    .data([], String);
    remove.enter().append(".matrix_canvas")
    ;
    
    remove.exit().remove();
    
}

function remove_matrix_marks(){
    remove = gene_model_display.selectAll(".matrix_mark")
    .data([], String);
    remove.enter().append(".matrix_mark")
    ;
    
    remove.exit().remove();
}

function forward_control_bar(){
    
    // alert("Max height bar: " + max_height_bar_genome_view);
    
    
    var control_bar_width = 200;
    var control_x_range = d3.scale.linear() // not sure what this is
    .domain([0, 18 - 1])
    .range([0, control_bar_width]); 
    
    
    // The bar we brush
    // This bar histogram needs to contain the bins.
    // each bin shows number of variants in this range
    
    
    // Have filtered variants - preserved, too! So have row entry
    // can bucket by range
    var num_bins = 18; // was 25
    var quality_range = forward_HQ_max;
    var interval = quality_range/num_bins;
    var interval_start = 0;
    var interval_end = interval_start+interval;
    var binned_counts_hist = [];
    var total_vars_possible = variants_with_exon_coords.length;
    for(i = 0; i < num_bins; i++){
        
        var filts = variants_with_exon_coords.filter(function(val) {
                                                     return ((val[34] <=interval_end)&&(val[34] >= interval_start));
                                                     });
        
        
        binned_counts_hist.push(filts.length);
        interval_start = interval_end;
        interval_end = interval_end + interval;
        
    }
    
    // Set the height to the largest count
    var max_count = Math.max.apply( Math, binned_counts_hist );
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(max_count)
    // .attr("text-anchor", "middle")
    .attr("x", 248)
    .attr("y", 20)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(0)
    // .attr("text-anchor", "middle")
    .attr("x", 248)
    .attr("y", 60)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    
    controlBar = controls_panel.selectAll("rect.itemsControl")
    .data(binned_counts_hist)
    .enter().append("svg:rect")
    .attr("class", "itemsControl")
    .attr("x", function(d, i) {return control_x_range(i) + 36;})
    .attr("y", function(d) {return (1 + 40 + 15 - (d/max_count)*40);})
    .attr("width",  control_bar_width/binned_counts_hist.length-2)
    // Bars come in upside down, can't flip them, so ...
    .attr("height", function(d) {return (d/max_count)*40}) // Note that changing this only changes the control
    // bar height
    .attr("fill", "#666");
    
    // The brush region we move
    // Need to pin down at one region and allow to fill
    controls_panel.append("g") // ? g What is this? This is the shaded brush itself
    .attr("class", "brush")
    .call(d3.svg.brush().x(d3.scale.linear().range([35, 35 + control_bar_width + control_bar_width/binned_counts_hist.length-2]))
          .extent([s_Forward_HQ[0],s_Forward_HQ[1]])
          .on("brush", brush_forward_HQ))
    .selectAll("rect")
    .attr("y", 16)
    .attr("height", 40);
    
}


function reverse_control_bar(){
    // Next Controller
    // exome
    
    var reverse_bar_width = 200;
    var forward_x_range = d3.scale.linear() // not sure what this is
    .domain([0, 18 - 1])
    .range([0, reverse_bar_width]); 
    
    
    // The bar we brush
    // This bar histogram needs to contain the bins.
    // each bin shows number of variants in this range
    
    
    // Have filtered variants - preserved, too! So have row entry
    // can bucket by range
    var num_bins = 18; // was 25
    var quality_range = reverse_HQ_max;
    var interval = quality_range/num_bins;
    var interval_start = 0;
    var interval_end = interval_start+interval;
    var binned_counts_hist = [];
    var total_vars_possible = variants_with_exon_coords.length;
    for(i = 0; i < num_bins; i++){
        var filts = variants_with_exon_coords.filter(function(val) {
                                                     return ((val[35] <=interval_end)&&(val[35] >= interval_start));
                                                     });
        
        
        binned_counts_hist.push(filts.length);
        interval_start = interval_end;
        interval_end = interval_end + interval;
        
    }
    
    // Set the height to the largest count
    var max_count = Math.max.apply( Math, binned_counts_hist );
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(max_count)
    // .attr("text-anchor", "middle")
    .attr("x", 497)
    .attr("y", 20)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(0)
    // .attr("text-anchor", "middle")
    .attr("x", 497)
    .attr("y", 60)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    
    rev_control_bar = controls_panel.selectAll("rect.itemsRevControl")
    .data(binned_counts_hist)
    .enter().append("svg:rect")
    .attr("class", "itemsRevControl")
    .attr("x", function(d, i) {return forward_x_range(i) + 250 + 36;})
    .attr("y", function(d) {return (1 + 40 + 15 - (d/max_count)*40);})
    .attr("width",  reverse_bar_width/binned_counts_hist.length-2)
    // Bars come in upside down, can't flip them, so ...
    .attr("height", function(d) {return (d/max_count)*40}) // Note that changing this only changes the control
    // bar height
    .attr("fill", "#666");
    
    // The brush region we move
    // Need to pin down at one region and allow to fill
    controls_panel.append("g") // ? g What is this? This is the shaded brush itself
    .attr("class", "brush")
    .call(d3.svg.brush().x(d3.scale.linear().range([250 + 35, 450 + 35 + reverse_bar_width/binned_counts_hist.length-2]))
          .extent([s_Reverse_HQ[0],s_Reverse_HQ[1]])
          .on("brush", brush_reverse_HQ))
    .selectAll("rect")
    .attr("y", 16)
    .attr("height", 40);
}


function ratio_control_bar(){
    
    // Need to create some ratios!
    
    // Calculate these up front!
    
    // index  32 33 34 35 Each is (34+35)/(32+ 33+ 34+ 35)
    
    // Putting the ratio value at the end for easy access
    
    
    
    var reverse_bar_width = 200;
    var forward_x_range = d3.scale.linear() // not sure what this is
    .domain([0, 18 - 1])
    .range([0, reverse_bar_width]); 
    
    
    // The bar we brush
    // This bar histogram needs to contain the bins.
    // each bin shows number of variants in this range
    
    
    // Have filtered variants - preserved, too! So have row entry
    // can bucket by range
    var num_bins = 18; // was 25
    var quality_range = ratio_HQ_max;
    var interval =quality_range/num_bins;
    var interval_start = 0;
    var interval_end = interval_start+interval;
    var binned_counts_hist = [];
    var total_vars_possible = variants_with_exon_coords.length;
    for(i = 0; i < num_bins; i++){
        var filts = variants_with_exon_coords.filter(function(val) {
                                                     return ((val[65] <=interval_end)&&(val[65] >= interval_start));
                                                     });
        
        
        binned_counts_hist.push(filts.length);
        interval_start = interval_end;
        interval_end = interval_end + interval;
    }
    
    // Set the height to the largest count
    var max_count = Math.max.apply( Math, binned_counts_hist );
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(max_count)
    // .attr("text-anchor", "middle")
    .attr("x", 748)
    .attr("y", 20)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    controls_panel.append("text")
    .attr("class", "hq_widget_range") 
    .text(0)
    // .attr("text-anchor", "middle")
    .attr("x", 748)
    .attr("y", 60)
    .attr("font-family", "sans-serif")
    .attr("font-size", "13px")
    .attr("fill", "#666");
    
    
    rev_control_bar = controls_panel.selectAll("rect.itemsRatioControl")
    .data(binned_counts_hist)
    .enter().append("svg:rect")
    .attr("class", "itemsRatioControl")
    .attr("x", function(d, i) {return forward_x_range(i) + 500 + 36;})
    .attr("y", function(d) {return (1 + 40 + 15 - (d/max_count)*40);})
    .attr("width",  reverse_bar_width/binned_counts_hist.length-2)
    // Bars come in upside down, can't flip them, so ...
    .attr("height", function(d) {return (d/max_count)*40}) // Note that changing this only changes the control
    // bar height
    .attr("fill", "#666");
    
    // The brush region we move
    // Need to pin down at one region and allow to fill
    controls_panel.append("g") // ? g What is this? This is the shaded brush itself
    .attr("class", "brush")
    .call(d3.svg.brush().x(d3.scale.linear().range([500 + 36, 700 + 36 + reverse_bar_width/binned_counts_hist.length-2]))
          .extent([s_Ratio_HQ[0],s_Ratio_HQ[1]])
          .on("brush", brush_ratio_HQ))
    .selectAll("rect")
    .attr("y", 16)
    .attr("height", 40);
    
}

// Callback function
function brush_ratio_HQ(){
    s_Ratio_HQ = d3.event.target.extent();
    // We might have to set globals        
    recurrence_histogram_bar(); // This could be where we set the bar
    
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    
    remove_matrix_marks();
    
    draw_matrix_marks();
    
}

// Callback function
function brush_forward_HQ() {
    
    s_Forward_HQ = d3.event.target.extent();
    // We might have to set globals        
    recurrence_histogram_bar(); // This could be where we set the bar
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    
    remove_matrix_marks();
    
    draw_matrix_marks();
}

// Callback function
function brush_reverse_HQ() {
    
    s_Reverse_HQ = d3.event.target.extent();
    
    recurrence_histogram_bar(); // This could be where we set the bar
    
    // Now change the variants onscreen
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    
    remove_matrix_marks();
    
    draw_matrix_marks();
}

// The large summary bar shows the culmination of all filters applied

// >>> Need recurrence percentage - need the number of patients to calculate this

function recurrence_histogram_bar() {
    
    
    // Wipe current slider values - they're old
    controls_panel.selectAll("text.hq_widget_values").remove();    
    
    var quality_range_reverse = reverse_HQ_max;
    var low_end_reverse = quality_range_reverse*Number(s_Reverse_HQ[0]);
    var high_end_reverse = quality_range_reverse*Number(s_Reverse_HQ[1]);
    
    
    controls_panel.append("text")
    .attr("class", "hq_widget_values") 
    .text("Range:  " + Math.floor(low_end_reverse) + " - " + Math.floor(high_end_reverse) + " reads")
    // .attr("text-anchor", "middle")
    .attr("x", 286)
    .attr("y", 72)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    
    // Reset
    filtered_variants = [];
    
    // The quality filters
    filtered_variants = variants_with_exon_coords.filter(function(val) {
                                                         return ((Number(val[35]) <=high_end_reverse)&&(Number(val[35]) >= low_end_reverse));
                                                         });
    
    var quality_range = forward_HQ_max; // For HQ could get dynamically
    var low_end = quality_range*Number(s_Forward_HQ[0]);
    var high_end = quality_range*Number(s_Forward_HQ[1]);
    
    
    controls_panel.append("text")
    .attr("class", "hq_widget_values") 
    .text("Range:  " + Math.floor(low_end) + " - " + Math.floor(high_end) + " reads")
    // .attr("text-anchor", "middle")
    .attr("x", 36)
    .attr("y", 72)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    
    filtered_variants = filtered_variants.filter(function(val) {
                                                 return ((Number(val[34]) <=high_end)&&(Number(val[34]) >= low_end));
                                                 });
    
    
    // x coord 536
    
    var quality_range_ratio = ratio_HQ_max; // For HQ could get dynamically
    var low_end_ratio = quality_range_ratio*Number(s_Ratio_HQ[0]);
    var high_end_ratio = quality_range_ratio*Number(s_Ratio_HQ[1]);
    
    
    controls_panel.append("text")
    .attr("class", "hq_widget_values") 
    .text("Range:  " + Math.floor(low_end_ratio) + " - " + Math.floor(high_end_ratio) + " %")
    .attr("x", 536)
    .attr("y", 72)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    
    filtered_variants = filtered_variants.filter(function(val) {
                                                 return ((Number(val[65]) <=high_end_ratio)&&(Number(val[65]) >= low_end_ratio));
                                                 });
    
    
    // Apply the check box filters
    
    
    if(nonsynonymous_filter_val == 0){
        
        filtered_variants = filtered_variants.filter(function(val) {
                                                     return (val[11].indexOf("NON_SYNONYMOUS_CODING(") == -1);
                                                     });
        
    }
    
    if(frame_shift_val == 0){
        
        filtered_variants = filtered_variants.filter(function(val) {
                                                     return (val[11].indexOf("FRAME_SHIFT(") == -1);
                                                     });
    }
    
    
    if(codon_deletion_filter_val == 0){
        
        filtered_variants = filtered_variants.filter(function(val) {
                                                     return (val[11].indexOf("CODON_DELETION(") == -1);
                                                     });
    }
    
    if(truncating_filter_val == 0){
        
        filtered_variants = filtered_variants.filter(function(val) {
                                                     return (val[11].indexOf("STOP_GAINED(") == -1);
                                                     });
    }
    
    // Have filtered variants - preserved, too! So have row entry
    // can bucket by range
    var num_bins = 200;
    
    var interval = Math.ceil(exon_model_length/num_bins);
    var interval_start = 0;
    var interval_end = 0+interval;
    var binned_counts_hist = [];
    var total_vars_possible = variants_with_exon_coords.length;
    for(i = 0; i < num_bins; i++){
        // Filtering on position
        var filts = filtered_variants.filter(function(val) {
                                             return ((Number(val[64]) <=interval_end)&&(Number(val[64]) >= interval_start));
                                             });
        
        
        binned_counts_hist.push(filts.length);
        interval_start = interval_end;
        interval_end = interval_end + interval;
        
    }
    
    var x_range = d3.scale.linear()
    .domain([0, binned_counts_hist.length - 1])
    .range([0, screen_width_gene_model]); 
    // wiping the bars currently there off the screen
    controls_panel.selectAll("rect.items").remove();
    // now add the new bars
    
    // scale to max value
    var max_count = Math.max.apply( Math, binned_counts_hist );
    
    if( filtered_variants.length !=0){
        
        // the data variable becomes "range"
        controls_panel.selectAll("rect.items")
        .data(binned_counts_hist)
        .enter().append("svg:rect")
        .attr("class", "items")
        .attr("x", function(d, i) {return x_range(i) + gene_display_x_offset;})
        // Change the height here //
        .attr("y", function(d) {return 1+ 240 - (d/number_patients)*90;})
        .attr("width",  screen_width_gene_model/200)
        .attr("fill", "#666")
        .attr("height", function(d) {return (d/number_patients)*90;});
    } else{
        
    }
    
    controls_panel.selectAll("text.total_num_variants_text").remove();
    total_num_variants = filtered_variants.length;
    
    controls_panel.append("text")
    .attr("class", "total_num_variants_text") 
    .text(total_num_variants)
    // .attr("text-anchor", "middle")
    .attr("x", 70)
    .attr("y", 198)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    // Calculate the total unique variants
    var list = [];
    
    for(i = 0; i < filtered_variants.length; i++){
        list.push(filtered_variants[i][5]);
    }
    var uniqueNames = [];
    $.each(list, function(i, el){
           if($.inArray(el, uniqueNames) === -1) uniqueNames.push(el);
           });
    
    list = [];
    var unique = uniqueNames.length;
    uniqueArray = [];
    controls_panel.selectAll("text.total_unique_variants_text").remove();
    controls_panel.append("text")
    .attr("class", "total_unique_variants_text") 
    .text(unique)
    // .attr("text-anchor", "middle")
    .attr("x", 70)
    .attr("y", 228)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    // Now change the variants onscreen
    
    
    
}
var total_num_variants;

var nonsynonymous_filter_val;

function nonsynonymous_filter(){
    recurrence_histogram_bar();
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    remove_matrix_marks();
    
    draw_matrix_marks();
}

function frame_shift_filter(){
    recurrence_histogram_bar();
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    remove_matrix_marks();
    
    draw_matrix_marks();
}

var codon_deletion_filter_val;

function codon_deletion_filter(){
    recurrence_histogram_bar();
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    remove_matrix_marks();
    
    draw_matrix_marks();
}

var truncating_filter_val;

function truncating_filter(){
    recurrence_histogram_bar();
    remove_current_variants();
    remove_protein_display();
    remove_exon_summary_view();
    draw_filtered_variants();
    drawProteinSummaryView();
    drawExonSummaryView();
    remove_matrix_marks();
    
    draw_matrix_marks();
}

function HQ_Forward_Control(){
    
}

function wipe_control_panel(){
    
    controls_panel.selectAll("rect.items").remove();
    controls_panel.selectAll("rect.itemsControl").remove();
    controls_panel.selectAll("rect.itemsRevControl").remove();
    controls_panel.selectAll("rect.itemsRatioControl").remove();
    controls_panel.selectAll("text.sample_text").remove();
    controls_panel.selectAll(".hq_widget_values").remove();
    controls_panel.selectAll(".hq_widget_range").remove(); 
    controls_panel.selectAll("g").remove();
}



// Coordinate conversion methods - these ensure that the 
// variant positions on each model can be comparable.

// The coordinate along a given transcript at which the variant
// will strike.
function computeAltTranscriptCoord(){
    
}

// The coordinate along the protein model at which the variant
// will strike.
function computeProteinModelCoord(){
    
}


function draw_filtered_variants(){
    // Variant Data
    for(i = 0; i < filtered_variants.length; i++){
        
        
        // Now need to translate to where the exon is on the screen
        var variant_pos_screen = (Number(filtered_variants[i][64])/exon_model_length)*screen_width_gene_model;
        
        
        gene_model_display
        .append("line")
        .attr("class", "variant_line_summ")
        .attr("x1", variant_pos_screen + gene_display_x_offset)
        .attr("y1", 2)
        .attr("x2", variant_pos_screen + gene_display_x_offset)
        .attr("y2", 355); // Bottom of screen is 355
        
        
        
    }
    
}

var filtered_samples_by_exon;

// The lines and text of the model
function drawPermanentStructures(){
}

// Global switch - don't invoke this method if switched on:
var alts_shown;

function drawAlternativeTrans(){
    
    var size_protein_display = 80;
    
    if(alts_shown == 0){
        // First lay down a canvas
        alts_shown = 10;
        var num_trans = 10;
        // The width might be okay (may also just need a complete redraw)
        // but should be able to calculate the length of the pane
        // on the fly.
        
        gene_model_display.append("rect")
        .attr("x", gene_display_x_offset)
        .attr("y", 70) // was 25 use this for exons
        .attr("width", 550)
        .attr("height", alts_shown*40 + size_protein_display) // Was 35 use this for exons
        .attr("stroke-width", 1.5)
        .attr("stroke", "#666")
        .attr("fill", "white")
        .on("mousedown", function(exon_id_number) 
            {   // Draw the alternative transcripts
            //d3.select(this).style("fill", "aliceblue");
            alts_shown = 0;
            //d3.select(this).style("opacity", 0);
            
            })
        .on("mouseup", function() 
            {
            // Remove all the alternative transcript exon blocks
            remove = gene_model_display.selectAll(".alt_trans")
            .data([], String);
            remove.enter().append("alt_trans")
            ;
            
            remove.exit().remove();
            
            // Remove all the alternative transcript splice sites
            remove = gene_model_display.selectAll(".alt_splices")
            .data([], String);
            remove.enter().append("alt_splices")
            ;
            
            remove.exit().remove();
            
            d3.select(this).style("fill", null);
            d3.select(this).remove();
            
            // Remove all the alternative transcript text names
            remove = gene_model_display.selectAll(".alt_trans_names")
            .data([], String);
            remove.enter().append("alt_trans_names")
            ;
            
            remove.exit().remove();
            
            d3.select(this).style("fill", null);
            d3.select(this).remove();
            
            });
        // Draw the alts
        // Call a function to do this:
        drawAlternativeTransModels();
        
    }else{
        
    }
}

function drawAlternativeTransModels(){
    // For now...
    numberAltTranscripts = 5;
    for(j=0;j < numberAltTranscripts;j++){
        // Put the id on each alt transcript:
        gene_model_display
        .append("text")
        .text("RUNX1")
        // .attr("text-anchor", "middle")
        .attr("x", 228)
        .attr("y", 85 + j*65)
        .attr("class","alt_trans_names");
        
        
        for(i = 0; i < numberOfExons; i++){
            var exon_id_number = i+1; // so we know which exon the user clicked on.
            gene_model_display
            .append("rect")
            .attr("class","alt_trans")
            .attr("x", gene_display_x_offset + exon_splicesite_screencoords[i])
            .attr("y", 90 + j*65) // was 25 use this for exons
            .attr("width", exonLengthsScreen[i])
            .attr("height", 35) // Was 35 use this for exons
            .attr("stroke-width", 1.5)
            .attr("stroke", "#E6550D")
            .attr("fill", "#FDAE6B")
            .on("mousedown", function(exon_id_number) 
                {   // Draw the alternative transcripts
                d3.select(this).style("fill", "aliceblue");
                
                })
            .on("mouseup", function() 
                {  
                d3.select(this).style("fill", null);  
                // Some sort of restore procedure... maybe fill it with null?
                // Try out "remove()" procedure?
                
                });
        }
        
        // Draw in the splice sites...
        for(i = 1; i < numberOfExons; i++){
            gene_model_display
            .append("rect")
            .attr("class","alt_splices")
            .attr("x", gene_display_x_offset + exon_splicesite_screencoords[i] - (length_of_splice_site*2/bp_length_model)*screen_width_gene_model)
            .attr("y", 95 + j*65) // was 25 use this for exons
            .attr("width", (length_of_splice_site*2/bp_length_model)*screen_width_gene_model)
            .attr("height", 25) // Was 35 use this for exons
            .attr("stroke-width", 1.5)
            .attr("stroke", "#E6550D")
            .attr("fill", "#E6550D");
        }
        
    }
    
    //alert("Fun");
}

// EXONS!

function remove_exon_summary_view(){
    remove = gene_model_display.selectAll(".exon_summ")
    .data([], String);
    remove.enter().append(".exon_summ")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll(".exon_text")
    .data([], String);
    remove.enter().append(".exon_text")
    ;
    
    remove.exit().remove();
}
var exon_view_y_offset = 110;

var track_stroke_width = 1.0;

var exon_spacer_offset;

//queried_gene_object


function drawExonSummaryView(){
    /*
     gene_model_display.append("text")
     .text("Transcript")
     .attr("x", 10)
     .attr("y", middle_pane_top_y_spacer + exon_view_y_offset -10)
     .attr("font-family", "sans-serif")
     .attr("font-size", "18px")
     .attr("fill", "#666");
     */
    exon_spacer_offset = exon_view_y_offset + 15;
    
    // Support for multiple transcripts
    
    /*
     gene_model_display
     .append("text")
     .text(queried_gene_object.gene_id)
     .attr("x", amino_acid_detail_x_offset + 5)
     .attr("y", middle_pane_top_y_spacer+exon_spacer_offset)
     .attr("font-family", "sans-serif")
     .attr("font-size", protein_annot_fontsize)
     .attr("fill", "#E6550D");
     */
    var accum = 0;
    var exon_lengths_total = 0;
    
    
    
    for(i = 0; i < queried_gene_object.exonLengths.length; i++){
        exon_lengths_total += queried_gene_object.exonLengths[i];
    }
    
    
    for(i = 0; i < queried_gene_object.exonLengths.length; i++){
        
        var screen_width_exon = screen_width_gene_model*(queried_gene_object.exonLengths[i]/exon_lengths_total);
        queried_gene_object.exon_lengths_total = exon_lengths_total;
        /*
         gene_model_display
         .append("rect")
         .attr("class", "exon_summ"  + String(j) + String(i))
         .attr("x", gene_display_x_offset + accum)
         .attr("y", middle_pane_top_y_spacer + exon_spacer_offset - 15) // was 25 use this for exons
         .attr("width", screen_width_exon)
         .attr("height", 15) // Was 35 use this for exons
         .attr("shape-rendering", "crispEdges" )
         .attr("stroke", "#E6550D")
         .attr("fill", "#FDAE6B")
         .attr("stroke-width", 0.75)
         .attr("fill-opacity", 0.75)
         .attr("stroke-opacity", 0.75);
         */
        accum += screen_width_exon;
    }
    
    
    
    
    // Put selectable, clear box on top to alow the user to choose one 
    
    var index = [j];
    
    var number_of_transcripts = g_all_gene_objects.length;
    
    gene_model_display
    .data(index)
    .append("rect")
    .attr("class", "transcript_number" + j)
    .attr("x", gene_display_x_offset )
    .attr("y", middle_pane_top_y_spacer + exon_spacer_offset - 15) // was 25 use this for exons
    .attr("width", screen_width_gene_model)
    .attr("height", 15) // Was 35 use this for exons
    .attr("stroke-width", 1.5)
    .attr("stroke", "#F0F0F0")
    .attr("fill", "#F0F0F0")
    .attr("fill-opacity", 0.0)
    .attr("stroke-opacity", 0.0);
    
    
}
// TOGO?
function drawExonSummary_All_Transcripts(){
    
    gene_model_display.append("text")
    .text("Transcript")
    .attr("x", 10)
    .attr("y", exon_view_y_offset -10)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    exon_spacer_offset = exon_view_y_offset + 15;
    
    // Support for multiple transcripts
    
    for(j = 0; j < g_all_gene_objects.length; j++){
        gene_model_display
        .append("text")
        .text(g_all_gene_objects[j].gene_id)
        .attr("x", amino_acid_detail_x_offset + 5)
        .attr("y", exon_spacer_offset)
        .attr("font-family", "sans-serif")
        .attr("font-size", protein_annot_fontsize)
        .attr("fill", "#E6550D");
        
        var accum = 0;
        var exon_lengths_total = 0;
        
        
        
        for(i = 0; i < g_all_gene_objects[j].exonLengths.length; i++){
            exon_lengths_total += g_all_gene_objects[j].exonLengths[i];
        }
        
        
        for(i = 0; i < g_all_gene_objects[j].exonLengths.length; i++){
            
            var screen_width_exon = screen_width_gene_model*(g_all_gene_objects[j].exonLengths[i]/exon_lengths_total);
            g_all_gene_objects[j].exon_lengths_total = exon_lengths_total;
            gene_model_display
            .append("rect")
            .attr("class", "exon_summ"  + String(j) + String(i))
            .attr("x", gene_display_x_offset + accum)
            .attr("y", exon_spacer_offset - 15) // was 25 use this for exons
            .attr("width", screen_width_exon)
            .attr("height", 15) // Was 35 use this for exons
            .attr("shape-rendering", "crispEdges" )
            .attr("stroke", "#E6550D")
            .attr("fill", "#FDAE6B")
            .attr("stroke-width", 0.75)
            .attr("fill-opacity", 0.25)
            .attr("stroke-opacity", 0.25);
            
            accum += screen_width_exon;
        }
        
        
        
        
        // Put selectable, clear box on top to alow the user to choose one 
        
        var index = [j];
        
        var number_of_transcripts = g_all_gene_objects.length;
        
        gene_model_display
        .data(index)
        .append("rect")
        .attr("class", "transcript_number" + j)
        .attr("x", gene_display_x_offset )
        .attr("y", exon_spacer_offset - 15) // was 25 use this for exons
        .attr("width", screen_width_gene_model)
        .attr("height", 15) // Was 35 use this for exons
        .attr("stroke-width", 1.5)
        .attr("stroke", "#F0F0F0")
        .attr("fill", "#F0F0F0")
        .attr("fill-opacity", 0.0)
        .attr("stroke-opacity", 0.0)
        .on("mouseover",function(d) {
            
            //  d3.select(this).style("stroke", "#000000");
            //  d3.select(this).style("stroke-opacity",1.0);
            
            // Select all associated variants
            //Hight light this transcript and all associated variants, fade out others
            
            for(i = 0; i < g_all_gene_objects[d].numberOfExons; i++){
            d3.select(".exon_summ"  +String(d) + String(i))
            .style("stroke", "#E6550D")
            .style("fill", "#FDAE6B")
            .style("stroke-width", 0.75)
            .style("fill-opacity", 0.75)
            .style("stroke-opacity", 1.0);
            }
            
            // Highlight all of their variants
            
            /*
             d3.select(".variant_triangle" + d).style("stroke", rollover_border_colouring).style("fill", rollover_border_colouring).style("stroke-opacity", 1.0);
             d3.select(".variant_line_summ" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
             
             
             d3.select(".mutation_type_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
             d3.select(".middle_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
             
             d3.select(".bottom_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
             */
            
            })
        .on("mouseout",function(d) {
            
            // d3.select(this).style("stroke", "#000000");
            // d3.select(this).style("stroke-opacity",0.0);
            
            for(i = 0; i < g_all_gene_objects[d].numberOfExons; i++){
            d3.select(".exon_summ" + String(d) + String(i))
            .style("stroke", "#E6550D")
            .style("fill", "#FDAE6B")
            .style("stroke-width", 0.75)
            .style("fill-opacity", 0.25)
            .style("stroke-opacity", 0.25);
            }
            
            
            
            
            /*
             // Fade all other out!
             for(i = 0; i < ){
             d3.select(".transcript_number" + d)
             .style("stroke", "#E6550D")
             .style("fill", "#FDAE6B")
             .style("stroke-width", 0.75)
             .style("fill-opacity", 0.25)
             .style("stroke-opacity", 0.25);
             }
             */
            /*
             //d3.select(".variant_box" + d).style("stroke", variant_objects[d].colour).style("stroke-opacity", 0.5);
             //d3.select(".bottom_box" + d).style("stroke", variant_objects[d].colour).style("stroke-opacity", 0.5);
             d3.select(".variant_triangle" + d).style("stroke", "#ddd").style("fill", "#ddd").style("stroke-opacity", 1.0);
             d3.select(".variant_line_summ" + d).style("stroke", "#ddd").style("fill", "#ddd").style("stroke-opacity", 1.0);
             
             
             d3.select(".bottom_box" + d).style("stroke", variant_objects[d].ref_aa_colour).style("stroke-opacity", 0.5);
             d3.select(".middle_box" + d).style("stroke", variant_objects[d].var_aa_colour).style("stroke-opacity", 0.5);
             
             
             
             d3.select(".mutation_type_box" + d).style("stroke",  "#C0B0FF").style("stroke-opacity", 0.5);
             
             */
            })
        .on("mousedown", function(d) 
            {
            // Leave highlighted permanently 
            
            
            })
        .on("mouseup", function(d) 
            {   
            
            })
        
        ;
        
        exon_spacer_offset += 20;
        
    }
}



// Protein Globals
var whole_chain;
var chain_name;

// Regions 
var protein_regions;
var region_coords;
var region_description;
//
// Sites
var site_coords;
var protein_sites;
//
// Modifications
var amino_acid_modifications;
var aa_modification_description;

function remove_protein_display(){
    remove = gene_model_display.selectAll(".protein")
    .data([], String);
    remove.enter().append(".protein")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll(".protein_desc_text")
    .data([], String);
    remove.enter().append(".protein_desc_text")
    ;
    
    remove.exit().remove();
}

var spacer;

var variant_objects;

var variant_indices;

var protein_annot_fontsize = "12px";

var previous_aa_coord;

var protein_annotations;
var protein_annot_colours;

function parse_aa_change_entry(p_aa_change){
    var ref_aas = "";
    var var_aas = "";
    
    var is_coord = 0;
    
    for(i = 0; i < p_aa_change.length; i++){
        if(p_aa_change[i] == '0' || p_aa_change[i] == '1' || p_aa_change[i] == '2' || p_aa_change[i] == '3'
           ||p_aa_change[i] == '4' || p_aa_change[i] == '5' || p_aa_change[i] == '6' || p_aa_change[i] == '7'
           || p_aa_change[i] == '8' || p_aa_change[i] == '9' ){
            
            is_coord = 1;
            
        }else{
            if(is_coord == 0){
                ref_aas += p_aa_change[i];
            }else{
                var_aas += p_aa_change[i];
            }
            
        }
    }
    
    var aa_change_object = new Object();
    aa_change_object.ref_aas = ref_aas;
    //alert(aa_change_object.ref_aas);
    aa_change_object.var_aas = var_aas;
    //alert(aa_change_object.var_aas);
    
    return aa_change_object;
}

var variant_line_stroke_colour = "#C2C2C2";


// Move globally - the aa colours
/*
 var electrically_charged = "#FF4019"; // Red
 var polar_uncharged = "#60E82E"; // Green
 var special_cases = "#9EABFF"; // blue/purple
 var hydrophobic_chain = "#F382FF"; // pink
 */

var electrically_charged = "#FF4019"; // Red
//var polar_uncharged = "#8BE602"; // Green
var polar_uncharged = "#8BE602"; //"#FFA85C"
//var special_cases = "#3333FF"; // blue/purple
var special_cases = "#7A7AFF"; 
var hydrophobic_chain = "#CCCCFF"; // pink




function drawProteinSummaryView(){
    variant_objects = [];
    
    // Background outlines, text etc...
    drawPermanentStructures(); 
    
    remove_variant_data();  
    
    gene_model_display.append("text")
    .text("Variants")
    .attr("x", 10)
    .attr("y", middle_pane_top_y_spacer + 15)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    
    // First, we need all the variants from every transcript - this will determine how much screen space we have.
    var all_variant_objects_coords = [];
    for(i = 0; i < g_all_gene_objects.length; i++){
        for(j = 0; j < g_all_gene_objects[i].variant_objects.length; j++){
            
            g_all_gene_objects[i].variant_objects[j].exome_screen_adjusted = g_all_gene_objects[i].variant_objects[j].exome_coord/g_all_gene_objects[i].exon_lengths_total;
            
            g_all_gene_objects[i].variant_objects[j].screen_id =  String(i) + String(j);
            all_variant_objects_coords.push(g_all_gene_objects[i].variant_objects[j]);
            
        }
    }
    
    var sorted_var_objects = all_variant_objects_coords.sort(function(a,b){
                                                             return a.exome_screen_adjusted - b.exome_screen_adjusted;
                                                             });
    
    var interm_sort_list = queried_gene_object.variant_objects.sort(function(a,b){
                                                                    return a.exome_coord - b.exome_coord;
                                                                    });
    
    queried_gene_object.variant_objects = interm_sort_list;
    
    variants_with_exon_coords = queried_gene_object.variants;
    // New coord sys. using aa coords directly from file.
    
    var aa_coords = [];
    var aa_changes = [];
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i++){
        
        
        aa_coords.push(queried_gene_object.variant_objects[i].aa_coord);
        aa_changes.push(queried_gene_object.variant_objects[i].aa_change);
        
        
        
    }
    for(j = 0; j < queried_gene_object.variant_objects.length; j++){
        
        var ref_aas = "";
        var var_aas = "";
        
        var is_coord = 0;
        var p_aa_change = queried_gene_object.variant_objects[j].aa_change;
        
        for(i = 0; i < p_aa_change.length; i++){
            if(p_aa_change[i] == '0' || p_aa_change[i] == '1' || p_aa_change[i] == '2' || p_aa_change[i] == '3'
               ||p_aa_change[i] == '4' || p_aa_change[i] == '5' || p_aa_change[i] == '6' || p_aa_change[i] == '7'
               || p_aa_change[i] == '8' || p_aa_change[i] == '9' ){
                
                is_coord = 1;
                
            }else{
                if(is_coord == 0){
                    ref_aas += p_aa_change[i];
                }else{
                    var_aas += p_aa_change[i];
                }
                
            }
        }
        
        
        queried_gene_object.variant_objects[j].ref_aas = ref_aas;
        
        queried_gene_object.variant_objects[j].var_aas = var_aas;
        
    }
    
    
    var max_end_index = 0;
    var max_index = 0;
    for(i = 0; i < whole_chains.length; i++){
        if(whole_chains[i][1] > max_end_index){
            max_end_index = whole_chains[i][1];
            max_index = i;
        }
    }
    
    var total_prot_model = max_end_index;
    var total_protein_model_screen = screen_width_protein_model*(total_prot_model/total_prot_model);
    
    
    // Grouping! 
    //  var scrn_width = screen_width_protein_model; // 675
    var scrn_width = screen_width_protein_model;
    
    var dist_with_scale_factor_sum = 0;
    
    var distances_between_aacoords = [];
    
    distances_between_aacoords.push(queried_gene_object.variant_objects[0].exome_coord); // the distance from the start of the model to the first variant.
    
    for(i = 0; i < queried_gene_object.variant_objects.length - 1; i++){
        var curr_dist = queried_gene_object.variant_objects[i+1].exome_coord - queried_gene_object.variant_objects[i].exome_coord; 
        distances_between_aacoords.push(curr_dist);
        
        dist_with_scale_factor_sum += curr_dist;
    }
    
    distances_between_aacoords.push(exon_model_length - queried_gene_object.variant_objects[queried_gene_object.variant_objects.length-1].exome_coord); // the distance from the last variant to the end of the model.
    
    var variant_block_width = 15;
    
    var num_variants = queried_gene_object.num_variants;
    
    var a_total = num_variants*variant_block_width;
    if(a_total > scrn_width){
        scrn_width = a_total;
        screen_width_protein_model = a_total;
    }else{
        
    }
    
    
    var reserved_for_variant_marks = num_variants*variant_block_width;
    
    var left_for_spacing = scrn_width - reserved_for_variant_marks;
    
    
    var scaled_spaced_between_variant_icons = [];
    
    for(i = 0; i < distances_between_aacoords.length; i++){
        var scaled_distance = (distances_between_aacoords[i]/exon_model_length)*left_for_spacing;
        scaled_spaced_between_variant_icons.push(scaled_distance);
    }
    
    var rollover_border_colouring = "#000000";
    
    variant_indices = [];
    var block_accum_spacing = scaled_spaced_between_variant_icons[0]; // Horiz spacing between variants
    
    
    // Figure out amino acid circle stacks formating
    // We need the biggest ref aa stack and the biggest var aa stack
    
    
    
    // find maxs
    var ref_max = 0;
    var var_max = 0;
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i++){
        
        if(ref_max < queried_gene_object.variant_objects[i].ref_aas.length){
            ref_max = queried_gene_object.variant_objects[i].ref_aas.length;
        }else{
            
        }
        
        if(var_max < queried_gene_object.variant_objects[i].var_aas.length){
            var_max = queried_gene_object.variant_objects[i].var_aas.length;
        }else{
            
        }
        
    }
    
    var aa_circle_size = variant_block_width;
    
    var var_aa_height = var_max*aa_circle_size; // the maximum height of the var aa stack.
    
    var ref_aa_height = ref_max*aa_circle_size; // the maximum height of the reference aa stack
    
    var pointer_arrow_pos = ref_aa_height;
    
    var pointer_arrow_height = 10;
    
    
    var transcript_and_protein_displacement = 0;
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i++){
        
        
        variant_indices.push(i);
        // Set their colour here to be reflective of effect
        // alert(variants_with_exon_coords[i][67]);
        var colour = "#E6550D";
        variant_object = new Object();
        variant_object.colour = colour;
        
        
        
        //}
        
        // The amino acid circle stacking
        var ref_aa_incrementer = middle_pane_top_y_spacer+35+7.5;
        var var_aa_incrementer = ref_aa_incrementer + pointer_arrow_pos + pointer_arrow_height;
        // Need to adjust so variants with fewer aa's start lower down the display
        var num_ref_aas = queried_gene_object.variant_objects[i].ref_aas.length;
        var spaces_left = ref_max - num_ref_aas;
        var add_to_incrementer = spaces_left*15;
        
        ref_aa_incrementer += add_to_incrementer;
        
        
        
        for(j = 0; j < queried_gene_object.variant_objects[i].ref_aas.length; j++){
            
            // Draw the vertical post to the mutation type:
            var num_aas = queried_gene_object.variant_objects[i].ref_aas.length;
            if((num_aas < ref_max) && (j==0)){
                
            }
            
            var aa_change_desc = aa_changes[i];
            
            var ref_aa = queried_gene_object.variant_objects[i].ref_aas[j];
            
            
            // Get the amino acid group colouring
            var ref_aa_colour = "#FFFFFF"; // default colour is white
            
            if(ref_aa == "R" || ref_aa == "H" || ref_aa == "K" || ref_aa == "D" || ref_aa == "E"){ 
                // then electrically-charged SC
                //ref_aa_colour = "#FF4019"; // Red
                ref_aa_colour = electrically_charged; // Red
            }else if(ref_aa == "S" || ref_aa == "T" || ref_aa == "N" || ref_aa == "Q" ){
                // then polar uncharged SC
                //ref_aa_colour = "#60E82E"; // Green
                ref_aa_colour = polar_uncharged; // Green
            }else if(ref_aa == "C" || ref_aa == "U" || ref_aa == "G" || ref_aa == "P" ){
                // Special Cases
                //ref_aa_colour = "#9EABFF"; // Blue/Purple
                ref_aa_colour = special_cases; // Blue/Purple
            }else if(ref_aa == "A" || ref_aa == "V" || ref_aa == "I" || ref_aa == "L" || ref_aa == "M"|| ref_aa == "F" || ref_aa == "Y" || ref_aa == "W"){
                // Hydrophobic Side Chain
                //ref_aa_colour = "#F382FF"; // Pink
                ref_aa_colour = hydrophobic_chain; // Yellow
            }else{
                ref_aa_colour = "#FFFFFF";
            }
            
            
            variant_object.ref_aa_colour = ref_aa_colour;
            
            variant_objects.push(variant_object);
            
            
            
            var stroke_circle = ref_aa_colour;
            if(ref_aa_colour ==  "#FFFFFF"){
                stroke_circle = "#000000" 
            }
            
            
            var index = [i]; // use this to index into global data structure
            gene_model_display
            .append("circle")
            .data(index)
            .attr("class", "bottom_box" + i)
            .attr("cx", gene_display_x_offset + block_accum_spacing + variant_block_width/2.0)
            .attr("cy",ref_aa_incrementer) // was 25 use this for exons
            .attr("r", variant_block_width/2.0)
            .attr("height", 15) // Was 25
            .attr("stroke-width", 1.5)
            .attr("stroke", stroke_circle)
            .attr("fill", ref_aa_colour)
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 1.0)
            .on("mousedown", function() 
                {   // Draw the alternative transcripts
                d3.select(this).style("fill", "aliceblue");
                
                })
            .on("mouseup", function() 
                {  
                d3.select(this).style("fill", null);  
                // Some sort of restore procedure... maybe fill it with null?
                // Try out "remove()" procedure?
                
                })
            .on("mouseover",function(d){
                
                })
            .on("mouseout",function(d){
                                });
            
            var text_aa_letter_adjust = 4;
            gene_model_display
            .append("text")
            .text(ref_aa)
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0)
            .attr("y", ref_aa_incrementer + text_aa_letter_adjust) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", protein_annot_fontsize)
            .attr("fill", "#000000");
            
            ref_aa_incrementer += 15;
        }
        
        // Draw the pointer/arrow triangle between the two stacks!
        //pointer_arrow_pos;
        var mid_point = "M " + (gene_display_x_offset + block_accum_spacing + 2) + " " + (ref_aa_incrementer-7.3) ;
        var triangle_point_2 = " L " + (gene_display_x_offset + block_accum_spacing + variant_block_width/2.0) + " " + (ref_aa_incrementer+10-7.5);
        var triangle_point_3 = " L " + (gene_display_x_offset + block_accum_spacing + variant_block_width - 2) + " " + (ref_aa_incrementer-7.3);
        
        var pointer_triangle = mid_point+triangle_point_2+triangle_point_3;
        
        
        gene_model_display
        .append("path")
        .attr("d", pointer_triangle + " z")
        .attr("stroke-width", 1.5)
        .attr("stroke", variant_line_stroke_colour)
        .attr("fill", variant_line_stroke_colour)
        .attr("fill-opacity", 0.75)
        .attr("stroke-opacity", 1.0);
        
        
        for(j = 0; j < queried_gene_object.variant_objects[i].var_aas.length; j++){
            
            var var_aa = queried_gene_object.variant_objects[i].var_aas[j];
            var var_aa_colour = "#FFFFFF"; // default colour is white
            
            if(var_aa == "R" || var_aa == "H" || var_aa == "K" || var_aa == "D" || var_aa == "E"){ 
                // then electrically-charged SC
                //var_aa_colour = "#FF4019"; // Red
                var_aa_colour = electrically_charged; // Red
            }else if(var_aa == "S" || var_aa == "T" || var_aa == "N" || var_aa == "Q" ){
                // then polar uncharged SC
                //var_aa_colour = "#60E82E"; // Green
                var_aa_colour = polar_uncharged; // Green
            }else if(var_aa == "C" || var_aa == "U" || var_aa == "G" || var_aa == "P" ){
                // Special Cases
                //var_aa_colour = "#9EABFF"; // Blue/Purple
                var_aa_colour = special_cases; // Blue/Purple
            }else if(var_aa == "A" || var_aa == "V" || var_aa == "I" || var_aa == "L" || var_aa == "M"|| var_aa == "F" || var_aa == "Y" || var_aa == "W"){
                // Hydrophobic Side Chain
                //var_aa_colour = "#F382FF"; // Pink
                var_aa_colour = hydrophobic_chain; // Yellow
            }else{
                var_aa_colour =  "#FFFFFF";
            }
        
            variant_object.var_aa_colour = var_aa_colour;
            var stroke_circle = var_aa_colour;
            if(var_aa_colour ==  "#FFFFFF"){
                stroke_circle = "#000000" 
            }
            
            gene_model_display
            .append("circle")
            .data(index)
            .attr("class", "middle_box" + i)
            .attr("cx", gene_display_x_offset + block_accum_spacing + variant_block_width/2.0)
            .attr("cy", var_aa_incrementer) // was 25 use this for exons
            .attr("r", variant_block_width/2.0)        
            .attr("stroke-width", 1.5)
            .attr("stroke", stroke_circle)
            .attr("fill", var_aa_colour)
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 1)
            
            .on("mouseover",function(d) {
               
                
                })
            .on("mouseout",function(d){
               
                });
            var text_aa_letter_adjust = 4;
            gene_model_display
            .append("text")
            .text(var_aa)
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0)
            .attr("y", var_aa_incrementer + text_aa_letter_adjust) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", protein_annot_fontsize)
            .attr("fill", "#000000");
            
            var_aa_incrementer += variant_block_width;
            
            
        }
        
        var exome_coord = queried_gene_object.variant_objects[i].exome_coord; // Make sure sorted!
        
        // Here we construct the path string!
        var triangle_aa_switch = (ref_aa_incrementer-7.3);
        
        var y_pos_triangle_path = triangle_aa_switch+var_max*variant_block_width + variant_block_width -4;
        
        var below_var_aa_spacer = y_pos_triangle_path + 40;
        
        var mid_point = "M " + (gene_display_x_offset + block_accum_spacing + 3) + " " + y_pos_triangle_path ;
        var triangle_point_2 = " L " + (gene_display_x_offset + screen_width_protein_model*(exome_coord/exon_model_length)) + " " + (below_var_aa_spacer);
        
        var triangle_point_3 = " L " + (gene_display_x_offset + block_accum_spacing + variant_block_width - 3) + " " + y_pos_triangle_path;
        
        var variant_triangle_path = mid_point+triangle_point_2+triangle_point_3;
        
        var index = [i]; 
        
        transcript_and_protein_displacement = below_var_aa_spacer;
        
        gene_model_display
        .append("path")
        .data(index)
        .attr("class", "variant_triangle" +i)
        .attr("d", variant_triangle_path + " z")
        .attr("stroke-width", 1.5)
        .attr("stroke", variant_line_stroke_colour)
        .attr("fill", variant_line_stroke_colour)
        .attr("fill-opacity", 0.75)
        .attr("stroke-opacity", 1.0)
        .on("mouseover",function(d) {
            d3.select(this).style("fill",  rollover_border_colouring);
            d3.select(this).style("stroke", rollover_border_colouring);
            d3.select(".variant_entry" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
            
            
            d3.select(".mutation_type_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
            // d3.select(".middle_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
            // d3.select(".variant_triangle" + d).style("stroke", rollover_border_colouring).style("fill", rollover_border_colouring).style("stroke-opacity", 1.0);
            // d3.select(".variant_line_summ" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
            // d3.select(".bottom_box" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0).style("cx", 10);
            
            })
        .on("mouseout",function(d) {
            d3.select(this).style("fill", variant_line_stroke_colour);
            d3.select(this).style("stroke", variant_line_stroke_colour);
            // d3.select(".bottom_box" + d).style("stroke","#FFFFFF").style("stroke-opacity", 0.5);
            // d3.select(".middle_box" + d).style("stroke","#FFFFFF").style("stroke-opacity", 0.5);
            
            // d3.select(".variant_line_summ" + d).style("stroke", "#ddd").style("fill", "#ddd").style("stroke-opacity", 1.0);                 
            
            d3.select(".mutation_type_box" + d).style("stroke",  variant_line_stroke_colour);
            d3.select(".variant_entry" + d).style("fill", "#F0F0F0");
            d3.select(".variant_entry" + d).style("stroke", "#F0F0F0");
            d3.select(".variant_entry" + d).style("stroke-opacity",0.50);
            
            });
        
        // Put a little post in to connect the bottom variant aa with the variant triangle below it!
        
        gene_model_display
        .append("rect")
        .attr("x", gene_display_x_offset + block_accum_spacing+3)
        .attr("y", var_aa_incrementer - variant_block_width/2.0) // was 25 use this for exons
        .attr("width", variant_block_width-3-3)
        .attr("height", y_pos_triangle_path-var_aa_incrementer+variant_block_width/2.0) // Was 25
        //.attr("height", 10)
        .attr("stroke-width", 1.5)
        .attr("stroke", variant_line_stroke_colour)
        .attr("fill", variant_line_stroke_colour)
        .attr("fill-opacity", 0.75)
        .attr("stroke-opacity", 1.0)
        
        // The top boxes.
        
        var index = [i]; // use this to index into global data structure
        
        
        // Harmless or Potential Cancer Causing: dbSNP and Cosmic visual encoding:
              
        var dbSNP_avail = false; 
        
        if(queried_gene_object.variant_objects[i].var_entry[dbSNP_129_index] != "" && queried_gene_object.variant_objects[i].var_entry[dbSNP_129_index] != "." && queried_gene_object.variant_objects[i].var_entry[dbSNP_129_index] != "ABSENT"){
            dbSNP_avail = true;
        }
        
        if(queried_gene_object.variant_objects[i].var_entry[dbSNP_135_index] != "" && queried_gene_object.variant_objects[i].var_entry[dbSNP_135_index] != "." && queried_gene_object.variant_objects[i].var_entry[dbSNP_135_index] != "ABSENT"){
            dbSNP_avail = true;
        }
        //alert(queried_gene_object.variant_objects[i].var_entry[dbSNP_137_index]);
        if(queried_gene_object.variant_objects[i].var_entry[dbSNP_137_index] != "" && queried_gene_object.variant_objects[i].var_entry[dbSNP_137_index] != "." && queried_gene_object.variant_objects[i].var_entry[dbSNP_137_index] != "ABSENT"){
            dbSNP_avail = true;
        }
        
        var cosmic_avail = false; 
        
        if(queried_gene_object.variant_objects[i].var_entry[cosmic_index] != "" && queried_gene_object.variant_objects[i].var_entry[cosmic_index] != "." && queried_gene_object.variant_objects[i].var_entry[cosmic_index] != "ABSENT"){
            cosmic_avail = true;
        }
        
        // Icons 
        
        var x_spacer_db_snp_cosmic = gene_display_x_offset + block_accum_spacing + 7.25;
        
        if(dbSNP_avail && cosmic_avail){
            //alert("1");
            gene_model_display
            .append("circle")
            .attr("cx", x_spacer_db_snp_cosmic)
            .attr("cy", middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width - 5) // was 25 use this for exons
            .attr("r", 3)        
            .attr("stroke-width", 1.5)
            .attr("stroke",  "#000000")
            .attr("fill",  "#000000")
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 0.50);
            
            gene_model_display
            .append("circle")
            .attr("cx", x_spacer_db_snp_cosmic)
            .attr("cy", middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width - 5 - 6) // was 25 use this for exons
            .attr("r", 3)        
            .attr("stroke-width", 1.5)
            .attr("stroke",  "#000000")
            .attr("fill",  "#FFFFFF")
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 0.50);
            
        }else if(dbSNP_avail){
           // alert("2")
            gene_model_display
            .append("circle")
            .attr("cx", x_spacer_db_snp_cosmic)
            .attr("cy", middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width - 5) // was 25 use this for exons
            .attr("r", 3)
            .attr("stroke-width", 1.5)
            .attr("stroke",  "#000000")
            .attr("fill",  "#FFFFFF")
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 0.50);
            
        }else if(cosmic_avail){
           // alert("cosmic");
            gene_model_display
            .append("circle")
            .attr("cx", x_spacer_db_snp_cosmic)
            .attr("cy", middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width - 5) // was 25 use this for exons
            .attr("r", 3)        
            .attr("stroke-width", 1.5)
            .attr("stroke",  "#000000")
            .attr("fill",  "#000000")
            .attr("fill-opacity", 0.50)
            .attr("stroke-opacity", 0.50);
            
        }
        
        gene_model_display
        .append("rect")
        .data(index)
        .attr("class", "mutation_type_box" + i)
        .attr("x", gene_display_x_offset + block_accum_spacing)
        .attr("y", middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width) // was 25 use this for exons
        .attr("width", variant_block_width)
        .attr("height", 15) // Was 25
        .attr("stroke-width", 2.5)
        .attr("stroke", variant_line_stroke_colour)
        .attr("fill", "#FFFFFF")
        .attr("fill-opacity", 0.50)
        .attr("stroke-opacity", 1.0)
        .on("mouseover",function(d) {
            
            d3.select(this).style("stroke", rollover_border_colouring);
            d3.select(this).style("stroke-opacity",1.0);
            
            d3.select(".variant_entry" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0);
            
            
            d3.select(".variant_triangle" + d).style("stroke", rollover_border_colouring).style("fill", rollover_border_colouring).style("stroke-opacity", 1.0);
            // d3.select(".variant_line_summ" + d).style("stroke", rollover_border_colouring).style("stroke-opacity", 1.0); 
            
            
            })
        .on("mouseout",function(d){
            
            d3.select(this).style("stroke",variant_line_stroke_colour);
            
            d3.select(".variant_triangle" + d).style("stroke",variant_line_stroke_colour).style("fill", variant_line_stroke_colour).style("stroke-opacity", 1.0);
            
            d3.select(".variant_entry" + d).style("fill", "#F0F0F0");
            d3.select(".variant_entry" + d).style("stroke", "#F0F0F0");
            d3.select(".variant_entry" + d).style("stroke-opacity",0.50);
            // d3.select(".variant_line_summ" + d).style("stroke", "#ddd").style("fill", "#ddd").style("stroke-opacity", 1.0); 
            });
        
        // Middle Variant AA Circle * Now the reference circles!
        
        
        
        
        // alert("Block Width = " + variant_block_width);
        var icon_symbol_y = middle_pane_top_y_spacer+20+ ref_max*variant_block_width - queried_gene_object.variant_objects[i].ref_aas.length*variant_block_width + 3;        
        var text_adjustment = 10;
        // Testing 
        if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("STOP_GAINED") != -1){
            gene_model_display
            .append("rect")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0)
            .attr("y", icon_symbol_y) // was 25 use this for exons
            .attr("width", variant_block_width/2.0)
            .attr("height", 15/2.0) // Was 25
            .attr("stroke-width", 1.5)
            .attr("stroke", "#000000")
            .attr("fill", "#000000")
            .attr("fill-opacity", 1.0)
            .attr("stroke-opacity", 1.0);
        }else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("FRAME_SHIFT(") != -1){
            
            gene_model_display
            .append("text")
            .text(">>")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0 - 2)
            .attr("y", icon_symbol_y + 8) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", protein_annot_fontsize)
            .attr("fill", "#000000");
        }else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("CODON_CHANGE_PLUS_CODON_DELETION(") != -1){
            
            var mid_point = "M " + (gene_display_x_offset + block_accum_spacing + 2 + 1) + " " + (icon_symbol_y+text_adjustment-7.3 -1) ;
            var triangle_point_2 = " L " + (gene_display_x_offset + block_accum_spacing + variant_block_width/2.0) + " " + (icon_symbol_y+text_adjustment+10-7.5 -4);
            var triangle_point_3 = " L " + (gene_display_x_offset + block_accum_spacing + variant_block_width - 2 - 1) + " " + (icon_symbol_y+text_adjustment-7.3-1);
            
            var pointer_triangle = mid_point+triangle_point_2+triangle_point_3;
            
            
            gene_model_display
            .append("path")
            .attr("d", pointer_triangle + " z")
            .attr("stroke-width", 1.5)
            .attr("stroke", "#000000")
            .attr("fill", "#000000")
            .attr("fill-opacity", 0.75)
            .attr("stroke-opacity", 1.0);
            
                
        }else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("CODON_DELETION(") != -1){
            
            gene_model_display
            .append("text")
            .text("-")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0 +1)
            .attr("y",icon_symbol_y + text_adjustment) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", "20px")
            .attr("fill", "#000000");
        }else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("CODON_INSERTION(") != -1){
            
            gene_model_display
            .append("text")
            .text("+")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0 -1)
            .attr("y", icon_symbol_y + text_adjustment) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", "18px")
            .attr("fill", "#000000");
        }
        else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("NON_SYNONYMOUS_CODING(") != -1){
            
            gene_model_display
            .append("rect")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0)
            .attr("y", icon_symbol_y) // was 25 use this for exons
            .attr("width", (variant_block_width/2.0)/2.0)
            .attr("height", 15/2.0) // Was 25
            .attr("stroke-width", 1.5)
            .attr("stroke", "#000000")
            .attr("fill", "#FFFFFF")
            .attr("fill-opacity", 1.0)
            .attr("stroke-opacity", 1.0);
            
            gene_model_display
            .append("rect")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0 + (variant_block_width/2.0)/2.0)
            .attr("y", icon_symbol_y) // was 25 use this for exons
            .attr("width", (variant_block_width/2.0)/2.0)
            .attr("height", 15/2.0) // Was 25
            .attr("stroke-width", 1.5)
            .attr("stroke", "#000000")
            .attr("fill", "#000000")
            .attr("fill-opacity", 1.0)
            .attr("stroke-opacity", 1.0);
            
        } else if(queried_gene_object.variant_objects[i].var_entry[variant_type_index].indexOf("SPLICE") != -1){
            
            gene_model_display
            .append("text")
            .text("*")
            .attr("x", gene_display_x_offset + block_accum_spacing + (variant_block_width/2.0)/2.0 + 1)
            .attr("y", icon_symbol_y+15) // was 25 use this for exons
            .attr("font-family", "sans-serif")
            .attr("font-size", "20px")
            .attr("fill", "#000000");
        }
        
        
        
        block_accum_spacing += variant_block_width+scaled_spaced_between_variant_icons[i+1];
        
    }
    
    
    // Labels!
    
    gene_model_display
    .append("text")
    .text("Mutation Type")
    .attr("x", amino_acid_detail_x_offset + 5)
    .attr("y", middle_pane_top_y_spacer + 20+15)
    .attr("font-family", "sans-serif")
    .attr("font-size", protein_annot_fontsize)
    .attr("fill", "#666");
    
    gene_model_display
    .append("text")
    .text("Reference A.A.s")
    .attr("x", amino_acid_detail_x_offset + 5)
    .attr("y", middle_pane_top_y_spacer + 35+15)
    .attr("font-family", "sans-serif")
    .attr("font-size", protein_annot_fontsize)
    .attr("fill", "#666");
    
    
    gene_model_display
    .append("text")
    .text("Variant A.A.s")
    .attr("x", amino_acid_detail_x_offset + 5)
    .attr("y", middle_pane_top_y_spacer + 50+15)
    .attr("font-family", "sans-serif")
    .attr("font-size", protein_annot_fontsize)
    .attr("fill", "#666");
    
    
    /// Variant stuff above! 
    
    /// Protein stuff below!
    
    // spacer = exon_spacer_offset; // this will specify the y coord at which each protein component will be 
    
    
    spacer= transcript_and_protein_displacement;
    
    gene_model_display.append("text")
    .text("Transcript")
    .attr("x", 10)
    .attr("y",spacer+5)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    exon_spacer_offset = exon_view_y_offset + 15;
    
    // Support for multiple transcripts
    
    spacer += 20
    gene_model_display
    .append("text")
    .text(queried_gene_object.gene_id)
    //.text("trans-anon")
    .attr("x", amino_acid_detail_x_offset + 5)
    .attr("y", spacer+12)
    .attr("font-family", "sans-serif")
    .attr("font-size", protein_annot_fontsize)
    .attr("fill", "#E6550D");
    
    var accum = 0;
    var exon_lengths_total = 0;
    
    
    
    for(i = 0; i < queried_gene_object.exonLengths.length; i++){
        exon_lengths_total += queried_gene_object.exonLengths[i];
    }
    
    var transcript_fill_colour = "#FFB87A";
    var transcript_stroke_colour = "#E6550D";
    
    var exon_y_pos = spacer;
    for(i = 0; i < queried_gene_object.exonLengths.length; i++){
        
        var screen_width_exon = screen_width_protein_model*(queried_gene_object.exonLengths[i]/exon_lengths_total);
        queried_gene_object.exon_lengths_total = exon_lengths_total;
        gene_model_display
        .append("rect")
        .attr("class", "exon_summ"  + String(j) + String(i))
        .attr("x", gene_display_x_offset + accum)
        .attr("y", spacer) // was 25 use this for exons
        .attr("width", screen_width_exon)
        .attr("height", 15) // Was 35 use this for exons
        .attr("shape-rendering", "crispEdges" )
        //.attr("stroke", "#E6550D")
        //.attr("fill", "#FDAE6B")
        .attr("stroke", transcript_stroke_colour)
        .attr("fill", transcript_fill_colour)
        .attr("stroke-width", 1.0)
        .attr("fill-opacity", 0.75)
        .attr("stroke-opacity", 1.0);
        
        accum += screen_width_exon;
    }
    
    spacer += 5;
    gene_model_display.append("text")
    .text("Protein")
    .attr("x", 10)
    .attr("y", middle_pane_top_y_spacer+spacer)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    spacer += 5;
    
    var aa_text_x_offset =  15;
    
    var fill_colours = ["#A1D99B","#9ECAE1", "#BCBDDC","#FFEDA0","#FFEDA0", "#FF5C61","#A5A8D6","#9ECAE1","#A1D99B", "#BCBDDC","#FFEDA0","#FFEDA0", "#FF5C61","#A5A8D6","#31A354","#3182BD","#756BB1","#FEB24C","#FF5CD9","#FF1F1F","#666885","#3182BD","#31A354","#756BB1","#FEB24C","#FF5CD9","#FF1F1F","#666885"];
    
    // Note that AA chain and Cleaved chain are neglected from this list because they are taken
    // care of below in the loop structure
    var prot_annot_labels = ["Signals","Domains", "Regions","Comp. Biases","Topo. Domains", "Transmem.","Zinc-Fingers","Active Sites","NP Binding","Metal Bind.","Bindings","Mod. Residue","Lipids", "Carbohyd.","Disuf."];
    
    // Molec Proc Class
    var molecule_processing = ["Signals", "A.A. Chain", "Cleaved Ch."];
    var molecule_processing_fill = "#A1D99B"; // green
    var molecule_processing_stroke = "#31A354"; // green
    
    // Regions class
    var regions_class = ["Domains", "Regions","Comp. Biases","Topo. Domains", "Transmem.","Zinc-Fingers"];
    var regions_class_fill = "#9ECAE1"; // blue
    var regions_class_stroke = "#3182BD"; // blue
    
    // Sites Class
    
    var sites_class = ["Active Sites","NP Binding","Metal Bind.","Bindings"];
    var sites_class_fill =  "#BCBDDC"; // Purps
    var sites_class_stroke = "#756BB1"; // Purps
    
    // AA Modifications Site
    
    var mods_class = ["Mod. Residue","Lipids", "Carbohyd.","Disuf."];
    var mods_class_fill =  "#FFBDBD"; // Pink
    var mods_class_stroke = "#FF726E"; // Pink
        
    var top_of_protein_annots = middle_pane_top_y_spacer+spacer;
    var stroke_colours = [];
    for(i = 0; i < protein_annotations.length; i++){
        // The inner loop - each annotation track (A.A. Chain, Domains, Regions, etc...)
        
        // Just do this by the value of the index: for instance, i = 0 is chain ... blah blah
        var annot_stroke = "";
        var annot_fill = "";
                
        if(i >=0 && i <1){
            annot_stroke = molecule_processing_stroke;
            annot_fill =  molecule_processing_fill;
        }else{
            annot_stroke = regions_class_stroke;
            annot_fill = regions_class_fill;
        }
        stroke_colours.push(annot_stroke);
        if(protein_annotations[i].length > 0){
            for(j = 0; j < protein_annotations[i].length; j++){
                
                
                
                var number_vars_striking = aa_coords.filter(function(val) {
                                                            return (protein_annotations[i][j][0] <= val && protein_annotations[i][j][1] >= val);
                                                            });
                
                
                // Compress to screen coords
                
                var seq_start = screen_width_protein_model*(protein_annotations[i][j][0]/total_prot_model);
                var seq_end = screen_width_protein_model*(protein_annotations[i][j][1]/total_prot_model);
                
                var width = seq_end - seq_start;
                if(width == 0){
                    width = 1;
                }
                if(number_vars_striking.length > 0){
                    var index = [[i,j]];
                    gene_model_display
                    .append("rect")
                    .data(index)
                    .attr("class", "protein")
                    .attr("x", gene_display_x_offset + seq_start)
                    .attr("y", middle_pane_top_y_spacer+spacer) // was 25 use this for exons
                    .attr("width", width)
                    .attr("height", 15) // Was 35 use this for exons
                    .attr("stroke-width", 0.5)
                    //.attr("stroke", function(d) {return stroke_colours[d[0]];})
                    //.attr("fill", function(d) {return fill_colours[d[0]];})
                    .attr("stroke", annot_stroke)
                    .attr("fill", annot_fill)
                    .attr("fill-opacity", 0.75)
                    .attr("stroke-opacity", 0.75)
                    .attr("shape-rendering", "crispEdges")
                    .on("mouseover",function(d) {
                        d3.select(this).style("stroke", rollover_border_colouring);
                        d3.select(this).style("stroke-width",1.5);
                        // Tool Tip
                        
                        gene_model_display
                        .append("rect")
                        .attr("class", "domain_tooltip")
                        .attr("x", Number(d3.select(this).attr("x")))
                        .attr("y", Number(d3.select(this).attr("y")) + 15) // was 25 use this for exons
                        .attr("width", prot_annot_info[d[0]][d[1]].length*8)
                        .attr("height", 20) // Was 35 use this for exons
                        .attr("stroke-width", 1.0)
                        .attr("stroke", "#FFFFFF")
                        .attr("fill", "#FFFFFF")
                        .attr("fill-opacity", 0.75)
                        .attr("stroke-opacity", 0.75)
                        
                        gene_model_display
                        .append("text")
                        .text(prot_annot_info[d[0]][d[1]])
                        .attr("class", "domain_tooltip_text")
                        .attr("x", Number(d3.select(this).attr("x")) + 3)
                        .attr("y", Number(d3.select(this).attr("y")) + 15 + 15) // was 25 use this for exons
                        .attr("font-family", "sans-serif")
                        .attr("font-size", "15px")
                        .attr("fill", "#666");
                        
                        
                        })
                    .on("mouseout",function(d) {
                        d3.select(this).style("stroke-width",0.5);
                        d3.select(this).style("stroke-opacity",0.75);
                        d3.select(this).style("stroke",stroke_colours[d[0]]);
                        // Remove domain tooltip
                        
                        d3.select(".domain_tooltip").remove();
                        d3.select(".domain_tooltip_text").remove();
                        
                        });
                }else{
                    var index = [[i,j]];
                    gene_model_display
                    .append("rect")
                    .data(index)
                    .attr("class", "protein")
                    .attr("x", gene_display_x_offset + seq_start)
                    .attr("y", middle_pane_top_y_spacer+spacer) // was 25 use this for exons
                    .attr("width", width)
                    .attr("height", 15) // Was 35 use this for exons
                    .attr("stroke-width", 0.5)
                    .attr("stroke", variant_line_stroke_colour)
                    .attr("fill", variant_line_stroke_colour)
                    .attr("fill-opacity", 0.3)
                    .attr("shape-rendering", "crispEdges")
                    .attr("stroke-opacity", 1.0)
                    
                    .on("mousedown", function() 
                        {   // Draw the alternative transcripts
                        d3.select(this).style("fill", "aliceblue");
                        
                        
                        })
                    .on("mouseup", function() 
                        {  
                        d3.select(this).style("fill", null);  
                        // Some sort of restore procedure... maybe fill it with null?
                        // Try out "remove()" procedure?
                        
                        })
                    .on("mouseover",function(d) {
                        d3.select(this).style("stroke", rollover_border_colouring);
                        d3.select(this).style("stroke-width",1.5);
                        // Tool Tip
                        
                        gene_model_display
                        .append("rect")
                        .attr("class", "domain_tooltip")
                        .attr("x", Number(d3.select(this).attr("x")))
                        .attr("y", Number(d3.select(this).attr("y")) + 15) // was 25 use this for exons
                        .attr("width", 150)
                        .attr("height", 20) // Was 35 use this for exons
                        .attr("stroke-width", 1.5)
                        .attr("stroke", "#FFFFFF")
                        .attr("fill", "#FFFFFF")
                        .attr("fill-opacity", 0.75)
                        .attr("stroke-opacity", 0.75)
                        
                        gene_model_display
                        .append("text")
                        .text(prot_annot_info[d[0]][d[1]])
                        .attr("class", "domain_tooltip_text")
                        .attr("x", Number(d3.select(this).attr("x")) + 3)
                        .attr("y", Number(d3.select(this).attr("y")) + 15 + 15) // was 25 use this for exons
                        .attr("font-family", "sans-serif")
                        .attr("font-size", "15px")
                        .attr("fill", "#666");
                        
                        
                        })
                    .on("mouseout",function() {
                        d3.select(this).style("stroke-width",0.5);
                        d3.select(this).style("stroke", variant_line_stroke_colour);
                        d3.select(this).style("stroke-opacity",1.0);
                        d3.select(".domain_tooltip").remove();
                        d3.select(".domain_tooltip_text").remove();
                        
                        });
                }
                if(i == 0){
                    if(j == 0){
                        spacer = spacer + 15;
                        gene_model_display
                        .append("text")
                        .attr("class","protein_desc_text")
                        .text("A.A. Chain")
                        .attr("x", aa_text_x_offset)
                        .attr("y", middle_pane_top_y_spacer+spacer)
                        .attr("font-family", "sans-serif")
                        .attr("font-size", protein_annot_fontsize)
                        .attr("fill", "#31A354");
                    }else{
                        spacer = spacer + 15;
                        gene_model_display
                        .append("text")
                        .attr("class","protein_desc_text")
                        .text("Cleaved Ch.")
                        .attr("x", aa_text_x_offset)
                        .attr("y", middle_pane_top_y_spacer+spacer)
                        .attr("font-family", "sans-serif")
                        .attr("font-size", protein_annot_fontsize)
                        .attr("fill", "#31A354");    
                    }
                }
            }
            
            if(i !=0){
                spacer = spacer + 15;
                //var label = prot_annot_labels[i-1];
                var label = prot_annot_labels[i-1];
                var another_index = [i];
                gene_model_display
                .append("text")
                .data(another_index)
                .attr("class","protein_desc_text")
                .text(label)
                .attr("x", aa_text_x_offset)
                .attr("y", middle_pane_top_y_spacer+spacer)
                .attr("font-family", "sans-serif")
                .attr("font-size", protein_annot_fontsize)
                .attr("fill", annot_stroke);
            }
        }
        
    }
    
    var bottom_protein_annots = middle_pane_top_y_spacer+spacer;
    // alert(queried_gene_object.sense);
    
    // Drawing the variant strikes on the display
    var triangle_bottom = y_pos_triangle_path + 40;
    
    for(i = 0; i < queried_gene_object.variant_objects.length; i ++){
        
        // The variant exome and aa coord line
        var exome_coord = queried_gene_object.variant_objects[i].exome_coord; 
        var aa_coord = queried_gene_object.variant_objects[i].aa_coord;
        var trans_strike_extent = 135;
        if(aa_coord == 0 || aa_coord > total_prot_model){
            trans_strike_extent = exon_y_pos+15;  // Bottom of the exon model
        }
        
        
        // Line that intersects the exon model
        gene_model_display
        .append("line")
        .attr("x1", (gene_display_x_offset + screen_width_protein_model*(exome_coord/exon_model_length)))
        .attr("y1", triangle_bottom)
        .attr("x2",  (gene_display_x_offset + screen_width_protein_model*(exome_coord/exon_model_length)))
        .attr("y2", exon_y_pos+15)
        .attr("stroke-width", 1.5)
        .attr("fill", "none")
        .attr("stroke", variant_line_stroke_colour)
        .attr("stroke-opacity", 1.0);
        
        
        if(aa_coord != 0 && aa_coord < total_prot_model){ // Skip if zero!
            
            // Middle line
            var index = [i];
            gene_model_display
            .data(index)
            .append("line")
            .attr("x1", (gene_display_x_offset + screen_width_protein_model*(exome_coord/exon_model_length)))
            .attr("y1", exon_y_pos+15)
            .attr("x2",  (gene_display_x_offset + screen_width_protein_model*(aa_coord/total_prot_model)))
            .attr("y2", top_of_protein_annots)
            .attr("stroke-width", 2)
            .attr("fill", "none")
            .attr("stroke", variant_line_stroke_colour)
            .attr("stroke-opacity", 1.0);
            
            // Line that intersects the protein model
            gene_model_display
            .append("line")
            .attr("x1", (gene_display_x_offset + screen_width_protein_model*(aa_coord/total_prot_model)))
            .attr("y1", top_of_protein_annots)
            .attr("x2",  (gene_display_x_offset + screen_width_protein_model*(aa_coord/total_prot_model)))
            .attr("y2", bottom_protein_annots)
            .attr("stroke-width", 1.5)
            .attr("fill", "none")
            .attr("stroke", variant_line_stroke_colour)
            .attr("stroke-opacity", 1.0);
        }
    }
    variant_table_y_spacer = bottom_protein_annots; 
}
var variant_table_y_spacer;
// To list the variant data - need a function to remove all these too!

/*
 This function draws the table at the bottom of the main gene view that contains
 all of the peripherally supporting variant attributes.
 
*/



function generate_variant_list() {
    // There's overlap in the fields - Cut-off fields?
    var rollover_border_colouring = "#000000";
    spacer = 30;
    
    // middle_pane_top_y_spacer = 20;
    
    wipe_variant_table_panel();
    
    gene_model_display
    .append("text")
    .text("Variant Data")
    .attr("x", 10)
    .attr("y", variant_table_y_spacer+spacer)
    .attr("font-family", "sans-serif")
    .attr("font-size", "18px")
    .attr("fill", "#666");
    
    spacer += 20;
    var entry_spacer = 0
    var offset = 70;
    var headings = ["Patient ID", "Chr. Coord.","Ref Base","Var Base", "dbSNP129", "dbSNP135","dbSNP137","COSMIC","A.A. Chng.", "Gene", "RefSeq ID"];    
    for(i = 0; i < headings.length; i++){
        gene_model_display
        .append("text")
        .text(headings[i])
        .attr("class", "variant_text")
        .attr("x", 20 + entry_spacer)
        .attr("y", variant_table_y_spacer+spacer + 12)
        .attr("font-family", "sans-serif")
        .attr("font-size", "12px")
        .attr("fill", "#666");
        entry_spacer += offset;
    }

    spacer += 17

    for(i = 0; i < queried_gene_object.variant_objects.length; i++){
        
        var within_entry_spacer = 0
        var spacer_offset = 70;
        
        var entry_coords = [0,1,2,3,4,5,6,7,9,10,11];
        
        for(j = 0; j <  entry_coords.length; j++){
                       
            gene_model_display
            .append("text")
            .text(queried_gene_object.variant_objects[i].var_entry[entry_coords[j]].substring(0, 9)) // may be better way...
            .attr("class", "variant_text")
            .attr("x", 20 + within_entry_spacer)
            .attr("y", variant_table_y_spacer+spacer + 12)
            .attr("font-family", "sans-serif")
            .attr("font-size", "12px")
            .attr("fill", "#666");
            within_entry_spacer += spacer_offset;

            
        }
        
        // For rollover clear boxes
        var index = [i];
        gene_model_display
        .data(index)
        .append("rect")
        .attr("class", "variant_entry" + i)
        .attr("x", 20 )
        .attr("y", variant_table_y_spacer+spacer) // was 25 use this for exons
        .attr("width", 770)
        .attr("height", 15) // Was 35 use this for exons
        .attr("stroke-width", 1.5)
        .attr("stroke", "#F0F0F0")
        .attr("fill", "#F0F0F0")
        .attr("fill-opacity", 0.50)
        .attr("stroke-opacity", 0.50)
        .on("mouseover",function(d) {
            d3.select(this).style("fill", "#ddd");
            d3.select(this).style("stroke", rollover_border_colouring);
            d3.select(this).style("stroke-opacity",1.0);
            
            d3.select(".variant_triangle" + d).style("stroke", rollover_border_colouring).style("fill", rollover_border_colouring).style("stroke-opacity", 1.0);
            d3.select(".mutation_type_box" + d).style("stroke", rollover_border_colouring);
        })
        .on("mouseout",function(d) {
            d3.select(this).style("fill", "#F0F0F0");
            d3.select(this).style("stroke", "#F0F0F0");
            d3.select(this).style("stroke-opacity",0.50);
            d3.select(".variant_triangle" + d).style("stroke",variant_line_stroke_colour).style("fill", variant_line_stroke_colour).style("stroke-opacity", 1.0);
            d3.select(".mutation_type_box" + d).style("stroke",  variant_line_stroke_colour);
        });
    
        spacer += 17; // A little space between each
    } 
}

var sorted_gene_objects;

function sort_alpha_order(){
    
    gene_objects = gene_objects.sort(function(a, b) {
                                     var textA = a.gene_name.toUpperCase();
                                     var textB = b.gene_name.toUpperCase();
                                     return (textA < textB) ? -1 : (textA > textB) ? 1 : 0;
                                     });
}

function sort_hotspot_score(){
    gene_objects = gene_objects.sort(function(a,b){
                                    return b.hotspot_score - a.hotspot_score;                                    });
    
}

function sort_variant_count(){
    
    gene_objects = gene_objects.sort(function(a,b){
                                     return b.num_variants - a.num_variants;
                                     });
}

function wipe_top_sort_list_panel(){
    
    sortable_list_display.selectAll("rect").remove();
    sortable_list_display.selectAll("text").remove();
}

function wipe_gene_selection_panel(){
    
    gene_selection_list.selectAll("rect").remove();
    gene_selection_list.selectAll("text").remove();
    
}

function wipe_variant_table_panel(){
    
    variant_table.selectAll("rect").remove();
    variant_table.selectAll("text").remove();
}

// Switches to determine which sorting option has been selected.
var file_sort = 0; // the orginal order... not sure how to recover this...
var alpha_sort = 0;
var hotspot_sort = 0;
var variantcount_sort = 0;
var num_truncs_sort = 0;


function gene_selection_display(){

    wipe_gene_selection_panel();
    wipe_top_sort_list_panel();
    
    var gene_selection_x_offset = 3;
    var gene_selection_y_offset = 3;
    
    var button_height = 20;
    var button_width = 100;
   
    gene_selection_y_offset += 5;
    
    gene_selection_list.append("text")
    .text("Sort Gene Selection List By:")
    .attr("x", gene_selection_x_offset + 30)
    .attr("y", gene_selection_y_offset + 4)
    .attr("font-family", "sans-serif")
    .attr("font-size", "15px")
    .attr("fill", "#666");
    
    
    gene_selection_y_offset += 2;
    
    gene_selection_list.append("text")
    .text("Alpha")
    .attr("x", gene_selection_x_offset + 35)
    .attr("y", gene_selection_y_offset + 22)
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "#666");
    
    gene_selection_list.append("text")
    .text("Cluster Score")
    .attr("x", gene_selection_x_offset + 35 + 48)
    .attr("y", gene_selection_y_offset + 22)
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "#666");
    
    gene_selection_list.append("text")
    .text("Variant Count")
    .attr("x", gene_selection_x_offset + 35 + 50 + 90)
    .attr("y", gene_selection_y_offset + 22)
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "#666");
    var alpha_fill = "#F0F0F0";
    var alpha_stroke = "#F0F0F0";
    var alpha_fill_opacity = 0.0;
    if(alpha_sort == 1){
        alpha_fill = "#ddd";
        alpha_fill_opacity = 0.45;
        alpha_stroke = "#000"
    }else{
        
    }
    
    gene_selection_list
    .append("rect")
    .attr("x", gene_selection_x_offset +30 )
    .attr("y",gene_selection_y_offset + 10) // was 25 use this for exons
    .attr("width", 50)
    .attr("height", button_height) // Was 35 use this for exons
    .attr("stroke-width", 1.5)
    .attr("stroke",alpha_stroke)
    .attr("fill", alpha_fill)
    .attr("fill-opacity", alpha_fill_opacity)
    .attr("stroke-opacity", 1.0)
    .on("mouseover",function() {
        d3.select(this).style("fill", "#ddd");
        d3.select(this).style("stroke", "#000000");
        d3.select(this).style("stroke-opacity",1.0);
        })
    .on("mouseout",function() {
        d3.select(this).style("fill", "#F0F0F0");
        d3.select(this).style("stroke", "#F0F0F0");
        })
    .on("mousedown", function() 
        {   
        d3.select(this).style("fill",  "#ddd").style("fill-opacity", 0.5);
        alpha_sort = 1;
        hotspot_sort = 0;
        variantcount_sort = 0;
        
        })
    .on("mouseup", function() 
        {  
        sort_alpha_order();
        getCurrGeneSelection(gene_objects[1].gene_id); 
        
        });
    
    var hot_fill = "#F0F0F0";
    var hot_stroke = "#F0F0F0";
    var hot_fill_opacity = 0.0;
    
    if(hotspot_sort == 1){
        hot_fill = "#ddd";
        hot_fill_opacity = 0.45;
        hot_stroke = "#000";
    }else{
        
    }
    
    gene_selection_list
    .append("rect")
    .attr("x", gene_selection_x_offset +30+50 )
    .attr("y",gene_selection_y_offset + 10) // was 25 use this for exons
    .attr("width", 90)
    .attr("height", button_height) // Was 35 use this for exons
    .attr("stroke-width", 1.5)
    .attr("stroke", hot_stroke)
    .attr("fill", hot_fill)
    .attr("fill-opacity", hot_fill_opacity)
    .attr("stroke-opacity", 1.0)
    .on("mouseover",function() {
        d3.select(this).style("fill", "#ddd");
        d3.select(this).style("stroke", "#000000");
        d3.select(this).style("stroke-opacity",1.0);
        })
    .on("mouseout",function() {
        d3.select(this).style("fill", "#F0F0F0");
        d3.select(this).style("stroke", "#F0F0F0");
        d3.select(this).style("stroke-opacity",0.50);
        })
    .on("mousedown", function() 
        {   // Draw the alternative transcripts
        d3.select(this).style("fill", "#ddd");
        
        
        })
    .on("mouseup", function() 
        {  
        d3.select(this).style("fill", "#F0F0F0");
        d3.select(this).style("stroke", "#F0F0F0");
        d3.select(this).style("stroke-opacity",0.50); 
        hotspot_sort = 1;
        alpha_sort = 0;
        variantcount_sort = 0;
        
        sort_hotspot_score();
        getCurrGeneSelection(gene_objects[1].gene_id); 
        // getCurrGeneSelection(gene_objects[d].gene_id); 
        
        });
    
    var count_fill = "#F0F0F0";
    var count_stroke = "#F0F0F0";
    var count_fill_opacity = 0.0;
    if(variantcount_sort == 1){
        count_fill = "#ddd";
        count_fill_opacity = 0.45;
        count_stroke = "#000"
    }else{
        
    }
    
    gene_selection_list
    .append("rect")
    .attr("x", gene_selection_x_offset +30+50+90 )
    .attr("y",gene_selection_y_offset + 10) // was 25 use this for exons
    .attr("width", 95)
    .attr("height", button_height) // Was 35 use this for exons
    .attr("stroke-width", 1.5)
    .attr("stroke", count_stroke)
    .attr("fill", count_fill)
    .attr("fill-opacity", count_fill_opacity)
    .attr("stroke-opacity", 1.0)
    .on("mouseover",function() {
        d3.select(this).style("fill", "#ddd");
        d3.select(this).style("stroke", "#000000");
        d3.select(this).style("stroke-opacity",1.0);
        })
    .on("mouseout",function() {
        d3.select(this).style("fill", "#F0F0F0");
        d3.select(this).style("stroke", "#F0F0F0");
        d3.select(this).style("stroke-opacity",0.50);
        })
    .on("mousedown", function() 
        {   // Draw the alternative transcripts
        d3.select(this).style("fill",  "#ddd");
        alpha_sort = 0;
        hotspot_sort = 0;
        variantcount_sort = 1;
        
        })
    .on("mouseup", function() 
        {  
        d3.select(this).style("fill", "#F0F0F0");
        d3.select(this).style("stroke", "#F0F0F0");
        d3.select(this).style("stroke-opacity",0.50); 
        sort_variant_count();
        getCurrGeneSelection(gene_objects[1].gene_id); 
        
        });
    
    
    var spacer = 20;
    
    var x_col_offset = 0;
    gene_selection_y_offset = 50;
    gene_selection_x_offset = 0;
    
    for(i = 0; i < gene_objects.length; i++){
        
        var index = [i];
        
        gene_selection_list.append("text")
        .text(gene_objects[i].gene_name + " (" + gene_objects[i].gene_id + ") ")
        .attr("class", "gene_selection_text" + i)
        .attr("x", gene_selection_x_offset + 35)
        .attr("y", spacer+gene_selection_y_offset + 4)
        .attr("font-family", "sans-serif")
        .attr("font-size", "12px")
        .attr("fill", "#666");
        
        gene_selection_list
        .data(index)
        .append("rect")
        .attr("class", "gene_selection_button" + i)
        .attr("x", gene_selection_x_offset +30 )
        .attr("y", spacer+gene_selection_y_offset -10) // was 25 use this for exons
        .attr("width", 155)
        .attr("height", button_height) // Was 35 use this for exons
        .attr("stroke-width", 1.5)
        .attr("stroke", "#F0F0F0")
        .attr("fill", "#F0F0F0")
        .attr("fill-opacity", 0.0)
        .attr("stroke-opacity", 1.0)
        .on("mouseover",function(d) {
            d3.select(this).style("fill", "#ddd");
            d3.select(this).style("stroke", "#000000");
            d3.select(this).style("stroke-opacity",1.0);
            })
        .on("mouseout",function(d) {
            d3.select(this).style("fill", "#F0F0F0");
            d3.select(this).style("stroke", "#F0F0F0");
            d3.select(this).style("stroke-opacity",0.50);
            })
        .on("mousedown", function() 
            {   // Draw the alternative transcripts
            d3.select(this).style("fill", "aliceblue");
            })
        .on("mouseup", function(d) 
            {  
            d3.select(this).style("fill", "#F0F0F0");
            d3.select(this).style("stroke", "#F0F0F0");
            d3.select(this).style("stroke-opacity",0.50); 
            
            getCurrGeneSelection(gene_objects[d].gene_id); 
            
            });
        spacer += 22;
    }
    
    
}

function remove_variant_data(){
    remove = gene_model_display.selectAll(".variant_entry")
    .data([], String);
    remove.enter().append(".variant_entry")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll(".variant_text")
    .data([], String);
    remove.enter().append(".variant_text")
    ;
    
    remove.exit().remove();
    
    remove = gene_model_display.selectAll(".variant_triangle")
    .data([], String);
    remove.enter().append(".variant_triangle")
    ;
    
    remove.exit().remove(); 
}

init_gene_model();
