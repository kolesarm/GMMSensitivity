
import os

def Main():
    # Instrument standard deviation table
    Parse('excluded_instrument_standard_deviation_matrix', 
          'excluded_instrument_std_dev')

    # Approximate sensitivity tables
    Parse('bootstrap_misspecification_table',     'bootstrap_table')
    
    # Sample sensitivity tables
    Parse('bootstrap_misspecification_table_sample',     
          'bootstrap_table_sample')

    # Global sensitivity table
    Parse('bootstrap_misspecification_table_global',     
          'bootstrap_table_global')
    
    # Global sensitivity table
    Parse('bootstrap_misspecification_table_global_sample',     
          'bootstrap_table_global_sample')
    
def Parse(intable, outtable):
    '''
    This function writes the tsv misspecification tables to new csv files
    '''
    infile  = '../output/' + intable  + '.tsv'
    outfile = '../output/' + outtable + '.csv'

    with open(infile, 'rU') as f:
        content = f.readlines()

    content = [x.strip(' \t').strip(' \n') for x in content]

    if intable in ["bootstrap_misspecification_table_global",
                   "bootstrap_misspecification_table_global_sample"]:
        content[len(content) - 1] = '      	'.join(content[len(content) - 1].split()[0:2])
        content[len(content) - 2] = '      	'.join(content[len(content) - 2].split()[0:2])

    with open(outfile, "wb") as f:
        for i in range(len(content)):
            if i == 0:
                f.write(content[i].strip())
            if i > 0:
                f.write(content[i])
            f.write('\n')

Main()
