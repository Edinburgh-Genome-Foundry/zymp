from zymp import (stacked_sites_array, plot_sequence_sites,
                  annotate_enzymes_sites, write_record)

enzymes_names = [
    'AccI', 'AclI', 'AflII', 'AflIII', 'AgeI', 'ApaLI', 'AseI',
    'AvaI', 'BamHI', 'BanII', 'BlnI', 'BmtI', 'BsmI', 'BssHII',
    'DdeI', 'DraI', 'Eco47III', 'EcoRI', 'EcoRV', 'HindII',
    'HindIII', 'HinfI', 'HpaI', 'KpnI', 'MfeI', 'MluI',
    'MspA1I', 'MunI', 'NaeI', 'NcoI', 'NdeI', 'NheI', 'NotI',
    'NsiI', 'NspI', 'PstI', 'PvuI', 'PvuII', 'SacI', 'SacII',
    'SalI', 'ScaI', 'SfaNI', 'SnaBI', 'SpeI', 'SphI', 'SspI',
    'StyI', 'VspI', 'XhoI', 'XmaI', 'ZraI'
]

forbidden_enzymes=['BsmBI', 'BsaI']

seq, sites_in_seq, leftover, success = stacked_sites_array(
        enzymes_names, forbidden_enzymes=forbidden_enzymes, tries=100)

print ("Sequence length:", len(seq),
       "\nRestriction sites:", len(sites_in_seq),
       "\nSites not included: ", leftover)
                  
# PLOT A SUMMARY
ax = plot_sequence_sites(seq, enzymes_names)
ax.figure.savefig("stacked_array.pdf", bbox_inches='tight')
                  
# WRITE THE SEQUENCE AND SITE ANNOTATIONS AS A RECORD
record = annotate_enzymes_sites(
    seq, enzymes_names, forbidden_enzymes=forbidden_enzymes)
write_record(record, 'stacked_site_array.gb')