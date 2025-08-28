from nupack import *
import random


def reverse_complement(sequence):
    complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in sequence[::-1]])


def generate_inhibition_strands():
    fixed_start = "GACUA" 
    fixed_end = "CUUUC" 
    non_conserved_seq = "GAGGCCGAAAGGCC"
    candidates = []

    for length in range(14, -1, -1):
        if length > 0:
            complement_part = reverse_complement(non_conserved_seq[:length])
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    mid_point = len(non_conserved_seq) // 2
    for length in range(14, -1, -1):
        if length > 0:
            start_idx = max(0, mid_point - length // 2)
            end_idx = min(len(non_conserved_seq), start_idx + length)
            complement_part = reverse_complement(non_conserved_seq[start_idx:end_idx])
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    for length in range(14, -1, -1):
        if length > 0:
            complement_part = reverse_complement(non_conserved_seq[-length:])
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    return candidates


def calculate_binding_energy(seq1, seq2, model):
    complex = Complex([Strand(seq1, name="strand1"), Strand(seq2, name="strand2")])
    result = complex_analysis(complexes=[complex], model=model, compute=['pfunc'])
    energy = result[complex].free_energy
    return energy


def generate_and_screen_irs(target_rna):
    model = Model(material='rna', celsius=31)
    recognition_strand = reverse_complement(target_rna)
    inhibition_candidates = generate_inhibition_strands()
    hibit_switch = "UCUCCUCUGGCGACCCUGAUGAGGCCGAAAGGCCGAAACGGUAUCGACCGUAGGUUGCCAGAACAGAGGAGAUAAAGAUGGUGAGCGGCUGGCGGCUGUUCAAGAAGAUUAGC"
    results = []
    for i, inhibition_strand in enumerate(inhibition_candidates):
        irs_sequence = inhibition_strand + recognition_strand
        dg1 = calculate_binding_energy(irs_sequence, target_rna, model)
        dg2 = calculate_binding_energy(irs_sequence, hibit_switch, model)
        displacement_score = abs(dg1) - abs(dg2)
        results.append({
            'irs_sequence': irs_sequence,
            'inhibition_strand': inhibition_strand,
            'recognition_strand': recognition_strand,
            'dg1': dg1,
            'dg2': dg2,
            'displacement_score': displacement_score
        })

    sorted_results = sorted(results, key=lambda x: x['displacement_score'], reverse=False)

    return sorted_results


def print_results(results, top_n=5):
    print("-" * 110)

    for i, result in enumerate(results[:top_n], 1):
        print(
            f"{i:<4} {result['irs_sequence']:<60} {result['dg1']:<15.2f} {result['dg2']:<15.2f} {result['displacement_score']:<15.2f}")

if __name__ == "__main__":
    target_rna_example = "GAUACUCCUAAUUAUGAUGUGCAGAAACACAUCAA"
    results = generate_and_screen_irs(target_rna_example)
    print_results(results, top_n=100)
    import csv

    with open('irs_screening_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['irs_sequence', 'inhibition_strand', 'recognition_strand',
                      'dg1', 'dg2', 'displacement_score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)



