import json

#table for each codon which other codon encodes same aa but is most distant
furthest_codon = dict([('GGT', 'GGC'), ('GGC', 'GGT'), ('TCG', 'AGC'), ('AGT', 'TCG'), ('AGC', 'AGT'), 
    ('ATG', 'ATG'), ('GCT', 'GCC'), ('GCC', 'GCT'), ('GCA', 'GCG'), ('GCG', 'GCA'), ('CGT', 'CGC'), ('CGC', 'CGT'),
    ('TGT', 'TGC'), ('TGC', 'TGT'), ('CTG', 'CTG'), ('CAA', 'CAG'), ('CAG', 'CAA'), ('CCT', 'CCG'), ('CCA', 'CCT'),
    ('CCG', 'CCA'), ('TAT', 'TAC'), ('TAC', 'TAT'), ('AAT', 'AAC'), ('AAC', 'AAT'), ('TTT', 'TTC'), ('TTC', 'TTT'),
    ('GTT', 'GTC'), ('GTC', 'GTT'), ('GTA', 'GTG'), ('GTG', 'GTA'), ('CAT', 'CAC'), ('CAC', 'CAT'), ('GAA', 'GAG'),
    ('GAG', 'GAA'), ('ACT', 'ACG'), ('ACC', 'ACT'), ('ACG', 'ACC'), ('ATT', 'ATC'), ('ATC', 'ATT'), ('TGG', 'TGG'),
    ('TAA', 'TGA'), ('TGA', 'TAA'), ('AAA', 'AAG'), ('AAG', 'AAA'), ('GAT', 'GAC'), ('GAC', 'GAT')])

json.dump( furthest_codon, open( "furthest_codon_pairs.json", 'w' ) )
