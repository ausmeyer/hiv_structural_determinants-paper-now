from Bio import SeqIO

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 
    'TGG':'W', '---':'-'
}

def translate_sequence(sequence):
  new_protein_sequence = ''
  new_dna_sequence = ''
  
  i = 0
  
  while i < len(sequence) - 1:
    codon = str(sequence[i]) + str(sequence[i+1]) + str(sequence[i+2])
    if codon in codontable.keys():
      new_dna_sequence += codon
      new_protein_sequence += codontable[codon]
    else:
      new_dna_sequence += '---'
      new_protein_sequence += '-'
    i += 3
    
  return(new_dna_sequence, new_protein_sequence)

def clean_records(records):
  keep_records = []
  base = ['A', 'C', 'T', 'G', 'a', 't', 'g', 'c', '-']
  for i, record in enumerate(records):
    keep = True
    for site in record:
      if site not in base:
        keep=False
        break
    if keep:
      keep_records.append(i)

  return(keep_records)
  
def main():
  records = list(SeqIO.parse(open('hiv1_protease.fasta', 'r'), 'fasta'))

  cleaned_records = clean_records(records)
  
  outfile_dna = open('hiv1_protease_clean_dna.fasta', 'w')
  outfile_protein = open('hiv1_protease_clean_protein.fasta', 'w')
  
  for record in cleaned_records:
    new_dna_sequence, new_protein_sequence = translate_sequence(str(records[record].seq))
    
    outfile_dna.write('>' + str(records[record].id) + '\n' + str(new_dna_sequence) + '\n')
    outfile_protein.write('>' + str(records[record].id) + '\n' + str(new_protein_sequence) + '\n')

#Run main program
if __name__ == '__main__':
   main()
