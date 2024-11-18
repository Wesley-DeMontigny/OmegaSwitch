from Bio import AlignIO
import sys

def fasta_to_nexus(input_fasta):
    alignment = AlignIO.read(input_fasta, "fasta")

    print("#NEXUS\n")
    print("BEGIN DATA;")
    print(f"DIMENSIONS ntax={len(alignment)} nchar={len(alignment[0])};")
    print("FORMAT datatype=dna gap=- missing=?;")
    print("MATRIX")

    for record in alignment:
        print(f"{record.id}\t\t\t\t\t\t\t\t\t{record.seq}")
    print(";")
    print("END;")


if __name__ == "__main__":
    input_fasta = sys.argv[1]
    fasta_to_nexus(input_fasta)
