import pysam
import argparse


"""
1) Correct aligned DddA BAM and replace likely deamination events with ambiguity codes (C|T: Y, G|A: R)
2) Create DA tag listing the moleular coordinates of deaminations
3) Add FD & LD tags for first and last deamination events in molecular coordinates
"""

# parse command line arguments
parser = argparse.ArgumentParser(description = "DddA BAM preprocessing",
    epilog = "")
parser.add_argument("-b", "--bam", required = True, metavar = '', help = "DddA aligned BAM to correct")
parser.add_argument("-c", "--cutoff", required = False, metavar = '', help = "Strand mut proportion cutoff")
args = parser.parse_args()

# identify fastq files in dir
bam_name = args.bam
if args.cutoff:
    cutoff = args.cutoff
else:
    cutoff = 0.90


def determine_da_strand_MD(read_obj, cutoff):
    # based on the proportion of C->T & G->A determine the strand acted upon by DddA
    # only counting single base substitutions
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    c = 0
    g = 0
    total = 0
    for pos in pair:
        if pos[0] == None or pos[1] == None: # indel, ignore
            pass
        else:
            qi = pos[0]
            ref_pos = pos[2].upper()
            if seq[qi] != ref_pos:
                total += 1
                change = ref_pos + seq[qi]
                if change == "CT":
                    c += 1
                elif change == 'GA':
                    g += 1
    if c+g == 0:
        return('none')
    elif c/(c+g) >= cutoff:
        return('CT')
    elif g/(c+g) >= cutoff:
        return('GA')
    else:
        return('undetermined')

def check_num_assigned(sam_obj, cutoff):
    none = 0
    und = 0
    ct = 0
    ga = 0
    for read in sam_obj.fetch():
        if read.is_secondary == False and read.is_supplementary == False:
            change = determine_da_strand_MD(read, cutoff)
            if change == 'none':
                none += 1
            elif change == 'undetermined':
                und += 1
            elif change == 'CT':
                ct += 1
            else:
                ga += 1
    return({'CT':ct, 'GA':ga, 'Undetermined':und, 'None':none})

def correct_read_MD(read_obj, strand):
    """ Identify single-base changes from the reference that are likely induced by DddA
     and correct the original sequence using ambiguity codes (C|T: Y, G|A: R) and output new DA-tag positions.
     Limit detection to the previously identified DddA strand info (either ct or ga).
     Output everything in FIBER coordinates, not reference!
    """
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    new_seq = ''
    amb_codes = {'CT':'Y', 'GA':'R'}
    deam_pos = [] # mol coordinates of likely base changes
    for pos in pair:
        if pos[0] == None: # deletion, ignore
            pass
        elif pos[1] == None: # insertion, use seq base
            qi = pos[0]
            new_seq += seq[qi]
        else:
            qi = pos[0]
            ref_pos = pos[2].upper()
            if seq[qi] != ref_pos:
                change = ref_pos + seq[qi]
                if change == strand:
                    new_seq += amb_codes[change] # Update seq with ambiguity codes
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
                else:
                    new_seq +=seq[qi]
            else:
                new_seq += seq[qi]
    return(new_seq, deam_pos)


# write corrected reads to new BAM
bam = pysam.AlignmentFile(bam_name, "rb")
new_bam = bam_name.replace('.bam','_corrected.bam')
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        strand = determine_da_strand_MD(read, cutoff)
        MD = read.get_tag('MD')
        if strand in ['CT','GA']:
            # WRITE NEW seq with added tags
            new_seq, deam_pos = correct_read_MD(read, strand)
            quals = read.query_qualities
            read.query_sequence = new_seq
            read.seq = new_seq
            read.query_qualities = quals
            read.set_tags([('DA', deam_pos), ('FD', deam_pos[0], "i"), ('LD', deam_pos[-1], "i"), ('ST', strand), ('MD', MD)])
        else:
            read.set_tags([('DA', [0]), ('FD', 0, "i"), ('LD', 0, "i"), ('ST', strand),('MD', MD)])
        corrected_bam.write(read)
bam.close()
corrected_bam.close()
