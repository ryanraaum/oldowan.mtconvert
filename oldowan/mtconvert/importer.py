import csv

def load_csv(file,
             header         = 1,     # number of rows to skip for header info
             hvr1           = False, # column number if present
             hvr1_covers    = [16024,16365],
             add16k         = True,  # add 16000 to every hvr1 site?
             hvr2           = False, # column number if present
             hvr2_covers    = [73,340],
             rflps          = False, # column number if present
             rflp_format    = False, # ? how to implement this ?
             sites          = False, # column number(s) if present
             sites_on_rCRS  = False, # matched entry or list for sites columns
             haplogroup     = False, # column number if present
             sample_id      = False, # column number if present
             haplotype_id   = False, # column number if present
             n              = False, # column number if present
             ):
    """Load mitochondrial haplotype definitions from csv file."""

    pass