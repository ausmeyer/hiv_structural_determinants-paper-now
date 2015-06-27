# open the file of new values (just 1 column of numbers, one for each alpha carbon)
inFile = open("distance.correlations", 'r')
 
# create the global, stored array
stored = []
 
# read the new B factors from file
for line in inFile.readlines(): stored.append( float(line) )
 
# close the input file
inFile.close()

max_b = max(stored)
min_b = min(stored)
 
# clear out the old B Factors
cmd.alter("4TVP_monomer_onlygp120 and n. CA", "b=0.0")
    
# update the B Factors with new properties
cmd.alter("4TVP_monomer_onlygp120 and n. CA", "b=stored.pop(0)")
 
# color the protein based on the new B Factors of the alpha carbons
cmd.spectrum("b", "rainbow_rev", "4TVP_monomer_onlygp120 and n. CA", minimum=min_b, maximum=max_b)
