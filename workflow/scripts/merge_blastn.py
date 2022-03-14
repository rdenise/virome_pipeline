import sys
import numpy as np
from numba import jit
import numpy.lib.recfunctions as rfn 
import os

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################
# Try to load it from file but not good now because can't access to the file
# if deployed worlflow by snakedeploy
##########################################################################

def iterrator_on_blast_hsp(blast_out):
    """Iterate over blast output. It considers that the output
    is in outfmt 6 and that all the hsp should be one after the
    other

    Args:
        blast_out (string): Path to a blast tabulated output

    Yields:
        pandas.DataFrame: A pandas Dataframe of the two proteins and their hsps
    """

    current_pair = ''
    list_line = []

    with open(blast_out) as r_file:

        for line in r_file:
            split_line = line.rstrip().split('\t')

            if (split_line[0], split_line[1]) == current_pair:
                list_line.append(split_line)
            elif current_pair == '':
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]
            else:
                yield list_line
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]

        yield list_line

##########################################################################

@jit(nopython=True)
def get_orientation(sstart, send):
    """Returns the orientation between two points.

    Args:
        sstart (numpy.array): Start of the sequence
        send (numpy.array): End of the sequence

    Returns:
        numpy.array: 1 or -1 depending on the orientation
    """
    
    return (send - sstart) // np.abs((send - sstart))

##########################################################################

@jit(nopython=True)
def update_values(old_values, sens):
    """Update sequence start/end with depending of the orientation.

    Args:
        old_values (numpy.array): Position in the sequence
        sens (numpy.array): 1 or -1 depending on the orientation

    Returns:
        np.array: Position in the sequence
    """

    return old_values * sens

##########################################################################

@jit(nopython=True)
def length_aln_on_sequence(start, end):
    """Return the length of the sequence.

    Args:
        start (numpy.array): Start of the sequence
        end (numpy.array): End of the sequence

    Returns:
        np.array: Length of the protein
    """
    return end - start + 1

##########################################################################

@jit(nopython=True)
def calculateIdentities(percIdentity, length):
    """Return the length of the sequence.

    Args:
        percIdentity (numpy.array): Percentage of identity from blast
        length (numpy.array): length of the alignment 

    Returns:
        np.array: number of identical amino acids in the alignment
    """    
    return np.floor(percIdentity * length / 100 + 0.5)

##########################################################################

@jit(nopython=True)
def calculate_fraction(delta, lgHSP, pid, pos):
    """Calculate the new score, identities and positives of the hsp2.

    Args:
        delta (int): Difference of size between old and new hsp2 length
        lgHSP (int): size of the hsp before trimming
        bitscore (float): Score of the alignment
        pid (int): Number of identical in the alignment 
        pos (int): Number of positives in the alignment

    Returns:
        float, int, int, int: The updated values of the score, identities, positives and length
    """

    # Calculate new score, id and positive
    # Calculation: initial_value * (franction of length that has been preserved)
    fraction = 1 - (delta / lgHSP)
    
    new_id = np.floor(pid * fraction)
    new_pos = np.floor(pos * fraction)
    
    # Calculate new length
    new_length = lgHSP - delta
    
    # Set expect value to -1 : this value should not be used after
    # having changed HSPs boundaries
    # new_evalue = -1
    
    return new_id, new_pos, new_length

##########################################################################

def remove_overlap_query(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the query.

    Args:
        hsp1 (np.array): Blast values of the given hsp (best)
        hsp2 (np.array): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    if hsp2[7] < hsp1[7]:
        new_qstart = hsp2[7]
        new_qend = hsp1[7] -1
        delta = hsp2[8] - new_qend
        new_sstart = hsp2[9]
        new_send = hsp2[10] - delta
        
    elif hsp2[8] > hsp1[8]:
        new_qstart = hsp1[8] + 1
        new_qend = hsp2[8]
        delta = new_qstart - hsp2[7] 
        new_sstart = hsp2[9] + delta
        new_send = hsp2[10]
        
    # lgHSP: 17, bitscore: 11, id: 12, pos:13
    new_id, new_pos, new_length = calculate_fraction(delta=delta, 
                                                    lgHSP=hsp2[17],
                                                    pid=hsp2[12],
                                                    pos=hsp2[13])
        
    return {17:new_length, 10:new_send,
            8:new_qend, 7:new_qstart,
            9:new_sstart, 13:new_pos, 12:new_id}


##########################################################################

def remove_overlap_subject(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the subject.

    Args:
        hsp1 (pandas.Series): Blast values of the given hsp (best)
        hsp2 (pandas.Series): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    if hsp2[9] < hsp1[9]:
        new_sstart = hsp2[9]
        new_send = hsp1[9] -1
        delta = hsp2[10] - new_send
        new_qstart = hsp2[7]
        new_qend = hsp2[8] - delta
        
    elif hsp2[10] > hsp1[10]:
        new_sstart = hsp1[10] + 1
        new_send = hsp2[10]
        delta = new_sstart - hsp2[9] 
        new_qstart = hsp2[7] + delta
        new_qend = hsp2[8]
        
    # lgHSP: 17, bitscore: 11, id: 12, pos:13
    new_id, new_pos, new_length = calculate_fraction(delta=delta, 
                                                    lgHSP=hsp2[17],
                                                    pid=hsp2[12],
                                                    pos=hsp2[13])
        
    return {17:new_length, 10:new_send,
            8:new_qend, 7:new_qstart, 
            9:new_sstart, 13:new_pos, 12:new_id}

##########################################################################

def checkHSPS(hsp1, hsp2, HSPMIN=100):
    """compare two HSPS in the blast output

    Args:
        hsp1 (np.array): Blast values of the given hsp (best)
        hsp2 (np.array): Blast values of the given hsp (questioning)
        HSPMIN (int, optional): Minumum length of the HSP. Defaults to 100.

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    dict_update = {18:0}
    
    # Check if the hsp2 is in a different orientation than hsp1: 'sens': 14
    if hsp1[14] != hsp2[14]:
        # print(f'orientation wrong: {hsp1[14]} != {hsp2[14]}')
        return dict_update
    
    # Check is hsp2 inside hsp1 for the query sequence: 'qstart': 6, 'qend': 7,
    if hsp1[7] >= hsp2[7] and hsp1[8] <= hsp2[8]:
        # print(f'hsp2 inside hsp1 for query: {hsp1[6]} >= {hsp2[6]} and {hsp1[7]} <= {hsp2[7]}')
        return dict_update
    
    # Check is hsp1 inside hsp2 for the query sequence: 'qstart': 6, 'qend': 7,
    elif hsp1[7] <= hsp2[7] and hsp1[8] >= hsp2[8]:
        # print(f'hsp1 inside hsp2 for query: {hsp1[6]} <= {hsp2[6]} and {hsp1[7]} >= {hsp2[7]}')
        return dict_update  

    # Check is hsp1 inside hsp2 for the subject sequence: 'sstart': 8, 'send': 9,
    elif hsp1[9] >= hsp2[9] and hsp1[10] <= hsp2[10]:
        # print(f'hsp2 inside hsp1 for subject: {hsp1[8]} >= {hsp2[8]} and {hsp1[9]} <= {hsp2[9]}')
        return dict_update
    
    # Check is hsp2 inside hsp1 for the subject sequence: 'sstart': 8, 'send': 9,
    elif hsp1[9] <= hsp2[9] and hsp1[10] >= hsp2[10]:
        # print(f'hsp1 inside hsp2 for subject: {hsp1[8]} <= {hsp2[8]} and {hsp1[9]} >= {hsp2[9]}')
        return dict_update 

    # reject HSPs that are in different orientation: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    # Query:    ---- A ---- B -----   A = HSP1 
    # Sbjct:    ---- B ---- A -----   B = HSP2

    if (hsp1[8] - hsp2[8]) * (hsp1[10] - hsp2[10]) < 0:
        # print(f'HSPs are in different orientation 1: ({hsp1[7]} - {hsp2[7]}) * ({hsp1[9]} - {hsp2[9]}) ===> {(hsp1[7] - hsp2[7]) * (hsp1[9] - hsp2[9])} < 0')
        return dict_update
    elif (hsp1[7] - hsp2[7]) * (hsp1[9] - hsp2[9]) < 0:
        # print(f'HSPs are in different orientation 2: ({hsp1[6]} - {hsp2[6]}) * ({hsp1[8]} - {hsp2[8]}) ===> {(hsp1[6] - hsp2[6]) * (hsp1[8] - hsp2[8])} < 0')
        return dict_update
    
    overlap_q = (hsp2[7] - hsp1[8]) * (hsp2[8] - hsp1[7])
    overlap_s = (hsp2[9] - hsp1[10]) * (hsp2[10] - hsp1[9])
    
    # Accept non-overlapping HSPs in correct orientation
    if overlap_q > 0 and overlap_s > 0 :
        # print(f'No overlap in query and subject: {overlap_q} > 0 and {overlap_s} > 0')
        dict_update[18] = 1
        return dict_update
    
    # Test if the query is overlaping
    elif overlap_q < 0:
        # print(f'Overlap in query: {overlap_q} > 0')
        dict_update = remove_overlap_query(hsp1=hsp1, 
                                           hsp2=hsp2)
        
        # update the hsp2 array with the new values
        for index_key in dict_update:
            hsp2[index_key] = dict_update[index_key]
            
        overlap_s = (hsp2[9] - hsp1[10]) * (hsp2[10] - hsp1[9])
    
    # Test if the subject is overlaping after possible update of an overlaping query
    if overlap_s < 0:
        # print(f'Overlap in subject: {overlap_s} > 0')
        dict_update = remove_overlap_subject(hsp1=hsp1, 
                                             hsp2=hsp2)
        
        # update the hsp2 array with the new values
        for index_key in dict_update:
            hsp2[index_key] = dict_update[index_key]
    
    # Filter out HSPs that are too short
    if hsp2[17] < HSPMIN:
        # print(f'HSP too short: {hsp2[17]} < {HSPMIN}')
        dict_update[18] = 0
        return dict_update
    
    # Set status to 1 for consistent HSPs
    dict_update[18] = 1
    return dict_update

##########################################################################

def prepare_df_hsps(list_hsps, blast_dtypes, HSPMIN=100):
    """Prepare a dataframe containing HSPs and update values.

    Args:
        list_hsps (list of strings): List of the lines containing the HSPs for the same pair of proteins
        blast_dtypes (dict): Type of the different columns in the blast outfmt 6 output
        HSPMIN (int, optional): [description]. Defaults to 100.

    Returns:
        pandas.DataFrame: Reduce DataFrame with only the significant hsps
    """    
    
    dtypes = np.dtype(blast_dtypes)
    
    df_hsps = np.array(list_hsps)
    df_hsps = rfn.unstructured_to_structured(df_hsps, dtype=dtypes)

    # Prepare the dtypes of the columns to add
    dtype2add = [
            ('id', np.int32),
            ('pos', np.int32),
            ('sens', np.int8),
            ('lg_aln_subject', np.int32),
            ('lg_aln_query', np.int32),
            ('lgHSP', np.int32),
            ('stat', np.int8),
           ]
    
    dtype2add = np.dtype(dtype2add)
    
    array2add = np.empty(df_hsps.shape, dtype=dtype2add)
        
    # Order HSPS by decreasing score
    # ind = np.argsort(df_hsps[:,11].astype(np.float64))[::-1]
    # ind = np.transpose([ind] * 18)
    # df_hsps = np.take_along_axis(df_hsps, ind, axis=0)
    
    # Sould add bitscore in output
    # ind = np.argsort(df_hsps, order='bitscore')[::-1]
    # df_hsps = df_hsps[ind]    
    
    # Calculate a more precise percentage of identity for each independent HSP: 'pident': 2, 'length': 3
    array2add['id'] = calculateIdentities(percIdentity=df_hsps['pident'], length=df_hsps['length'])
    array2add['pos'] = array2add['id']
    
    # Get the sens of the HSP on the subject: 'sstart': 8, 'send': 9
    array2add['sens'] = get_orientation(sstart=df_hsps['sstart'], send=df_hsps['send'])
    
    # Update the values of the position in case the sens is changed for the subject (as done in silix)
    df_hsps['sstart'] = update_values(old_values=df_hsps['sstart'], sens=array2add['sens'])
    df_hsps['send'] = update_values(old_values=df_hsps['send'], sens=array2add['sens'])
    
    # Calculate the length of the alignment on the subject and query: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    array2add['lg_aln_subject'] = length_aln_on_sequence(start=df_hsps['sstart'], end=df_hsps['send'])
    array2add['lg_aln_query'] = length_aln_on_sequence(start=df_hsps['qstart'], end=df_hsps['qend'])
    
    # Calculate length of the HSP: lg_aln_subject: 15, lg_aln_query:16 in silix it is max but I don't know why
    array2add['lgHSP'] = np.min([array2add['lg_aln_subject'], array2add['lg_aln_query']], axis=0)
    
    # Put Stats value to know which HSPS to conserve
    array2add['stat'] = (np.array([1] + [-1]*(df_hsps.shape[0] - 1)))
    
    # Concat df_hsps and array2add
    df_hsps = rfn.merge_arrays([df_hsps, array2add], flatten = True, usemask = False)
    
    for hsp1_index in range(df_hsps.shape[0] - 1):
        for hsp2_index in range(hsp1_index + 1, df_hsps.shape[0]):
            if df_hsps['stat'][hsp1_index] and df_hsps['stat'][hsp2_index]:
                values2update = checkHSPS(hsp1=df_hsps[hsp1_index], 
                                          hsp2=df_hsps[hsp2_index], 
                                          HSPMIN=HSPMIN)
                
                # print(hsp1_index, hsp2_index, values2update)
                for value in values2update:
                    df_hsps[hsp2_index][value] = values2update[value]
    
    return df_hsps[df_hsps['stat'] == 1]

##########################################################################

@jit(nopython=True)
def calculate_coverage(length_total, length_query, length_subject, option_cov='mean'):
    """Calculate the coverage of a sequence.

    Args:
        length_total (int): Length of the HSPs
        length_query (int): Length of the query
        length_subject (int): Length of the subject
        option_cov (str, optional): Silix option for calculating the coverage. Defaults to 'mean'.

    Returns:
        float: The percentage of coverage of the HSPs 
    """
    if option_cov == 'mean':
        cov = length_total / ((length_query + length_subject)/ 2.)
    elif option_cov == 'subject':
        cov = length_total / length_subject
    elif option_cov == 'query':
        cov = length_total / length_query
    elif option_cov == 'shortest':
        cov = length_total / min(length_query, length_subject)
    elif option_cov == 'longest':
        cov = length_total / max(length_query, length_subject)
    
    return cov

##########################################################################

@jit(nopython=True)
def calculate_percid(pos_total, length_query, length_subject, length_total, option_pid='mean'):
    """Calculates the percentage of identity of the total positives of the sequence.

    Args:
        pos_total (int): Number of positives in the alignment
        length_query (int): Length of the query
        length_subject (int): Length of the subject
        length_total (int): Length of the HSPs
        option_pid (str, optional): Silix option for the percentage of identitities calculation. Defaults to 'mean'.

    Returns:
        float: Percentage of identities in the alignment
    """
    if option_pid == 'mean':
        pid = pos_total / ((length_query + length_subject)/2.)
    elif option_pid == 'subject':
        pid = pos_total / length_subject
    elif option_pid == 'query':
        pid = pos_total / length_query
    elif option_pid == 'shortest':
        pid = pos_total / min(length_query, length_subject)
    elif option_pid == 'longest':
        pid = pos_total / max(length_query, length_subject)
    elif option_pid == 'HSP':
        pid = pos_total / length_total
        
    return pid

##########################################################################

def summarize_hits(df_hsps, option_cov='mean', option_pid='mean'):
    """Calculate summary statistics for a HSP

    Args:
        df_hsps (pandas.DataFrame): Sub dataframe with only the significative hsps
        length_query (int): Length of the query
        length_subject (int): Length of the subject
        option_cov (str, optional): Silix option for the coverage calculation. Defaults to 'mean'.
        option_pid (str, optional): Silix option for the percentage of identitities calculation. Defaults to 'mean'.

    Returns:
        int, float, float, float, float: Return the needed values to filter the hits
    """
    
    # Sort HSPs by query start position
    df_hsps = np.sort(df_hsps, order='qstart')
    
    # N-term delta between query and subject
    delta_lg = np.abs(df_hsps['qstart'] - df_hsps['sstart'])
    
    # C-term delta between query and subject and add to length difference between aligned sequences
    delta_lg += np.abs((df_hsps['qlen'][0] - df_hsps['qend'][-1]) - (df_hsps['slen'][0] - df_hsps['send']))
    
    # Internal gaps 
    query_starts = df_hsps['qstart'][1:]
    query_ends = df_hsps['qend'][:-1]
    subject_starts = df_hsps['sstart'][1:]
    subject_ends = df_hsps['send'][:-1]
    
    delta_lg += np.sum(np.abs((query_starts - query_ends) - (subject_starts - subject_ends)))

    # Calculate total length, score, Identities and Positives :
    length_total = np.sum(df_hsps['lgHSP'])
    pos_total = np.sum(df_hsps['pos'])
    # id_total = np.sum(df_hsps.id.values)
    
    # Cumulative length of HSPs / shortest seq. length = % coverage
    frac_HSPshlen = calculate_coverage(length_total=length_total, 
                                       length_query=df_hsps['qlen'][0], 
                                       length_subject=df_hsps['slen'][0], 
                                       option_cov=option_cov)
    
    # Number of positives / shortest seq. length = percentage of identity
    pospercentsh = calculate_percid(pos_total=pos_total, 
                                    length_query=df_hsps['qlen'][0], 
                                    length_subject=df_hsps['slen'][0], 
                                    length_total=length_total, 
                                    option_pid=option_pid)
    
    evalue = max(df_hsps['evalue'])

    return delta_lg, frac_HSPshlen, pospercentsh, evalue, df_hsps['stitle'][0]

##########################################################################

def summarize_hit_only(split_line, option_cov='mean', option_pid='mean'):
    """Summarize information of a single hit.

    Args:
        split_line (list of string): Split line of the blast output
        length_query (int): Length of the query
        length_subject (int): Length of the subject
        option_cov (str, optional): Silix option for the coverage calculation. Defaults to 'mean'.
        option_pid (str, optional): Silix option for the percentage of identitities calculation. Defaults to 'mean'.

    Returns:
        float, float, float, float: Return the needed values to filter the hits
    """

    # Get all the identical in the sequence: pident: 2, length: 3
    perc_id = float(split_line[2])
    aln_length = int(split_line[3])
    pos = np.floor(perc_id * aln_length / 100 + 0.5)

    qstart = int(split_line[7])
    qend = int(split_line[8])
    sstart = int(split_line[9])
    send = int(split_line[10])

    # Calculation as if multiple hit
    lg_aln_subject = length_aln_on_sequence(start=sstart, end=send)
    lg_aln_query = length_aln_on_sequence(start=qstart, end=qend)

    lg_total = max(lg_aln_query, lg_aln_subject)

    perc_id = calculate_percid(pos_total=pos,
                               length_query=int(split_line[4]),
                               length_subject=int(split_line[5]),
                               length_total=lg_total, 
                               option_pid=option_pid)

    cov = calculate_coverage(length_total=lg_total,
                             length_query=int(split_line[4]),
                             length_subject=int(split_line[5]), 
                             option_cov=option_cov)

    # Get the rest: evalue: 6
    evalue = float(split_line[6])

    return perc_id, cov, evalue, split_line[11]


##########################################################################
##########################################################################

blastn_files = snakemake.input.all_out

evalue = float(snakemake.wildcards.evalue)
coverage = float(snakemake.wildcards.coverage)

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'qstart', 'qend',
               'sstart', 'send', 'stitle']


# Get the types of te columns for multiple HSPs dataframe
blast_dtypes = [('qseqid','S100'),
                ('sseqid','S100'),
                ('pident',np.float64),
                ('length',np.int32),
                ('qlen',np.int32),
                ('slen',np.int32),
                ('evalue',np.float64),
                ('qstart',np.int32),
                ('qend',np.int32),
                ('sstart',np.int32),
                ('send',np.int32),
                ('stitle','S200'),
            ]


with open(snakemake.output.tsv, "w") as w_file:
    header = "\t".join(['protein1', 'protein2', 'pident', 'evalue', 'coverage', 'database', 'stitle'])
    w_file.write(f"{header}\n")

    num_files = len(blastn_files)

    for blast_index in range(num_files):
        blast_file = blastn_files[blast_index]

        database_file = blast_file.split('.')[-4]
        database = os.path.join(snakemake.params.dict_databases[database_file]['path'])

        for sub_blast in iterrator_on_blast_hsp(blast_out=blast_file) :
            # Get the number of hsps
            num_HSPs = len(sub_blast)

            qseqid = sub_blast[0][0]
            sseqid = sub_blast[0][1]

            if num_HSPs == 1:
                pident_blast, coverage_blast, evalue_blast, description = summarize_hit_only(split_line = sub_blast[0], 
                                                                                option_cov = snakemake.config['default_blast_option']['cov_min'],
                                                                                option_pid = snakemake.config['default_blast_option']['pid_min'])
            else:
                df_hsps = prepare_df_hsps(list_hsps = sub_blast,
                                        blast_dtypes = blast_dtypes, 
                                        HSPMIN = snakemake.params.minimum_length)

                if df_hsps.shape[0] == 1:
                    pident_blast, coverage_blast, evalue_blast, descrition = summarize_hit_only(split_line = df_hsps[0], 
                                                                                                option_cov = snakemake.config['default_blast_option']['cov_min'],
                                                                                                 option_pid = snakemake.config['default_blast_option']['pid_min'])
                else:
                    delta_lg, coverage_blast, pident_blast, evalue_blast, description = summarize_hits(df_hsps = df_hsps, 
                                                                                                        option_cov = snakemake.config['default_blast_option']['cov_min'], 
                                                                                                      option_pid = snakemake.config['default_blast_option']['pid_min'])


            if evalue_blast <= evalue and coverage >= coverage_blast :
                line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\t{database}\t{description}\n'
                w_file.write(line2write)


###########################################################
###########################################################

