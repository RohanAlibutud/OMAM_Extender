"""
   ___       _   _           __  __       __  __ 
  / _ \ _ __| |_| |__   ___ |  \/  | __ _|  \/  |
 | | | | '__| __| '_ \ / _ \| |\/| |/ _` | |\/| |
 | |_| | |  | |_| | | | (_) | |  | | (_| | |  | |
  \___/|_|   \__|_| |_|\___/|_|  |_|\__,_|_|  |_|
                                                 
  _____      _                 _                 
 | ____|_  _| |_ ___ _ __   __| | ___ _ __       
 |  _| \ \/ / __/ _ \ '_ \ / _` |/ _ \ '__|      
 | |___ >  <| ||  __/ | | | (_| |  __/ |         
 |_____/_/\_\\__\___|_| |_|\__,_|\___|_|         
 
 Updated 2025.01.09
                                                 
"""

#||---------------||
#||   LIBRARIES   ||
#||---------------||
import argparse
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


#||-------------||
#||   METHODS   ||
#||-------------||

# ARGUMENT PARSER
def parse_arguments():
    parser = argparse.ArgumentParser(description = "Extend a set of alignments with available species data from NCBI")
    parser.add_argument("--spp", required = True, help = "Species to be added")
    parser.add_argument("--omam", required = True, help = "Folder containing OrthoMaM alignments")
    parser.add_argument("--out", required = True, help = "Output folder for extended alignments")
    
    args = parser.parse_args()
    
    # convert to boolean
    args.extraction_on = args.extraction_on.lower() in ['true', '1', 'yes']
    
    return args

# METHOD TO BLAST A SEQUENCE AND RETRIEVE RESULTS
def blast_sequence(sequence_record):
    print("Running BLAST search...")
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_record.format("fasta"))
    blast_records = NCBIXML.read(result_handle)
    return blast_records

# METHOD TO ADD NEW SPECIES TO ALIGNMENT USING CONSENSUS SEQUENCE
def add_new_species(consensus_seq, records, spp_add):
    for species in spp_add:
        print(f"Retrieving sequence for {species} using consensus sequence...")
        blast_result = blast_sequence(consensus_seq)
        
        if not blast_result.alignments:
            print(f"No significant BLAST alignments found for {species}. Skipping.")
            continue
        
        for alignment in blast_result.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.001 and len(hsp.sbjct) >= 0.7 * len(consensus_seq):
                    sequence_str = hsp.sbjct.replace("-", "")
                    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
                    
                    if sequence_str and all(char in valid_amino_acids for char in sequence_str):
                        new_seq = Seq(sequence_str)
                        new_record = SeqRecord(
                            new_seq, 
                            id=species,
                            description=f"New sequence added from BLAST hit: {alignment.hit_def}"
                        )
                        records.append(new_record)
                        print(f"Added {species} to alignment with sequence:\n{new_record.seq}")
                    else:
                        print(f"Invalid sequence for {species}, skipping this hit.")
                    break
            else:
                continue
            break
        else:
            print(f"No high-quality alignments found for {species}. Skipping.")

    return records

# METHOD TO CALCULATE CONSENSUS SEQUENCE
def calculate_consensus(records):
    alignment = AlignIO.MultipleSeqAlignment(records)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    return consensus

# METHOD TO PROCESS AN INDIVIDUAL ALIGNMENT
def process_alignment(fasta_id, fasta_name, in_folder, out_folder, spp_add):
    print(f"Processing alignment: {fasta_id}")
    
    # Open alignment file and convert it to 2-line FASTA format
    records = list(SeqIO.parse(f"{in_folder}/{fasta_name}", "fasta"))
    
    # Calculate consensus sequence
    consensus_seq = calculate_consensus(records)
    
    # Convert consensus sequence to SeqRecord for BLAST formatting
    consensus_record = SeqRecord(consensus_seq, id="Consensus_Seq", description="Consensus sequence for BLAST")
    
    # Add new species to alignment
    updated_records = add_new_species(consensus_record, records, spp_add)
    
    # Save alignment in 2-line FASTA format
    out_file = f"{out_folder}/{fasta_name}"
    with open(out_file, "w") as output_handle:
        for record in updated_records:
            record.seq = Seq(str(record.seq).replace("\n", ""))  # Ensure each sequence is single-line
            SeqIO.write(record, output_handle, "fasta")
    
    print(f"Alignment saved to {out_file}")

#||--------------||
#||   COMMANDS   ||
#||--------------||
if __name__ == "__main__":
    args = parse_arguments()
    
    # Create output folder if it doesn't exist
    os.makedirs(args.out, exist_ok=True)
    
    # iterate over files in folder
    for file in os.listdir(args.omam):
        if file == "paths.txt":
            continue
        else:
            fasta_id = file.split("_")[1]
            process_alignment(
                fasta_id,
                file,
                args.omam,
                args.out,
                args.spp
                )
        
