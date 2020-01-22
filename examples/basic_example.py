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

forbidden_enzymes = ['BsmBI', 'BsaI']

sequence, enzymes_in_sequence, enzymes_not_in_sequence = stacked_sites_array(
    enzymes_names, forbidden_enzymes=forbidden_enzymes, tries=100
)

print (
    "Sequence length:", len(sequence),
    "\nRestriction sites:", len(enzymes_in_sequence),
    "\nSites not included: ", enzymes_not_in_sequence
)
                  
# PLOT A SUMMARY
ax = plot_sequence_sites(sequence, enzymes_names)
ax.figure.savefig("stacked_array.pdf", bbox_inches='tight')
                  
# WRITE THE SEQUENCE AND SITE ANNOTATIONS AS A RECORD
record = annotate_enzymes_sites(
    sequence, enzymes_names, forbidden_enzymes=forbidden_enzymes)
write_record(record, 'stacked_site_array.gb')