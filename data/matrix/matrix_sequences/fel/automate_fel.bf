/* Automates FEL with the following input. Written 10/1/2014 AGM. */

BASEDIR = "/home/austin/Desktop/hiv_structural_determinants/data/matrix/matrix_sequences/fel/";
datafile="hiv1_matrix_clean_dna.fasta";
output="run.log";
treefile="hiv1_matrix_clean_dna.tree";
sites="sites.dat";

inputRedirect = {};
inputRedirect["01"]="Universal";             //Genetic code
inputRedirect["02"]="New Analysis";          //New analysis
inputRedirect["03"]=BASEDIR+datafile;        //Fasta file, full path
inputRedirect["04"]="Default";               //Use HKY85 and MG94xHKY85
inputRedirect["05"]=BASEDIR+treefile;        //Tree
inputRedirect["06"]=BASEDIR+output;          //Output
inputRedirect["07"]="Estimate dN/dS only";   //Only estimate w
inputRedirect["08"]="One rate FEL";          //Estimate one rate
inputRedirect["09"]="0.1";                   //Significance level
inputRedirect["10"]=BASEDIR+sites;           //sitewise output

ExecuteAFile ("QuickSelectionDetection.bf", inputRedirect);
