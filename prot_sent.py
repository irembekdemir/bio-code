import json
import re

# (mRNA -> Amino Acid)
CODON_TABLE = {
    'AUG': 'Met (Start)', 'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

class BioinformaticsEngine:
    """Engine class performing nucleotide conversions and protein synthesis simulations."""

    COMPLEMENT_DNA = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    COMPLEMENT_RNA = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    TRANSCRIPTION_MAP = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}

    @staticmethod
    def validate_dna(dna_seq: str) -> tuple[bool, str]:
        """Validates the DNA sequence and detects invalid characters."""
        clean_seq = dna_seq.upper().strip()
        invalid_chars = set(re.sub(r'[ATCG]', '', clean_seq))
        if invalid_chars:
            return False, f"Invalid nucleotides detected: {', '.join(invalid_chars)}"
        return True, clean_seq

    @classmethod
    def replicate_dna(cls, dna_seq: str) -> dict:
        """DNA replication: Complementary strand (3'->5') and Reverse Complement (5'->3')."""
        is_valid, msg = cls.validate_dna(dna_seq)
        if not is_valid:
            return {"error": msg}

        complement = "".join(cls.COMPLEMENT_DNA[base] for base in msg)
        reverse_complement = complement[::-1]
        
        return {
            "original_5to3": msg,
            "complement_3to5": complement,
            "reverse_complement_5to3": reverse_complement
        }

    @classmethod
    def transcribe(cls, dna_coding_strand: str) -> dict:
        """mRNA synthesis (Transcription) from the DNA template strand."""
        is_valid, msg = cls.validate_dna(dna_coding_strand)
        if not is_valid:
            return {"error": msg}

        mrna = "".join(cls.TRANSCRIPTION_MAP[base] for base in msg)
        return {"mrna_sequence": mrna}

    @classmethod
    def translate_mrna(cls, mrna_seq: str) -> dict:
        """Translates mRNA into a polypeptide chain and tRNA anticodons at the ribosome."""
        codons = [mrna_seq[i:i+3] for i in range(0, len(mrna_seq), 3)]
        
        translation_results = []
        protein_chain = []
        
        for idx, codon in enumerate(codons, start=1):
            if len(codon) < 3:
                translation_results.append({
                    "codon_index": idx,
                    "mrna_codon": codon,
                    "trna_anticodon": "N/A",
                    "amino_acid": "Incomplete Codon (Partial)"
                })
                continue
                
            try:
                anticodon = "".join(cls.COMPLEMENT_RNA[base] for base in codon)
            except KeyError:
                anticodon = "Undefined Base"

            amino_acid = CODON_TABLE.get(codon, "Unknown Codon")
            protein_chain.append(amino_acid)

            translation_results.append({
                "codon_index": idx,
                "mrna_codon": codon,
                "trna_anticodon": anticodon,
                "amino_acid": amino_acid
            })

        return {
            "codons_analyzed": translation_results,
            "polypeptide_chain": " - ".join(protein_chain)
        }

    @classmethod
    def process_pipeline(cls, sample_id: str, dna_seq: str) -> dict:
        """Executes the full bioinformatics workflow."""
        is_valid, check_msg = cls.validate_dna(dna_seq)
        
        if not is_valid:
            return {
                "sample_id": sample_id,
                "status": "ERROR",
                "message": check_msg,
                "input_sequence": dna_seq
            }

        repl_data = cls.replicate_dna(check_msg)
        transcribe_data = cls.transcribe(check_msg)
        translate_data = cls.translate_mrna(transcribe_data["mrna_sequence"])

        return {
            "sample_id": sample_id,
            "status": "SUCCESS",
            "input_dna_5to3": check_msg,
            "replication": repl_data,
            "transcription_mrna": transcribe_data["mrna_sequence"],
            "ribosome_translation": translate_data
        }


def load_fasta(file_path: str) -> dict:
    """Reads the sample FASTA input file and returns a dictionary of sequences."""
    sequences = {}
    current_id = None
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line
                
    return sequences


if __name__ == "__main__":
  
    fasta_file = "sample_dna.fasta"
    data = load_fasta(fasta_file)

    pipeline_results = []
    for sample_id, dna_seq in data.items():
        res = BioinformaticsEngine.process_pipeline(sample_id, dna_seq)
        pipeline_results.append(res)

    print(json.dumps(pipeline_results, indent=2, ensure_ascii=False))
