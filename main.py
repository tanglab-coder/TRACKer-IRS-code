from nupack import *
import random


def reverse_complement(sequence):
    """生成RNA序列的反向互补序列"""
    complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in sequence[::-1]])


def generate_inhibition_strands():
    """按照特定规则生成候选抑制链序列"""
    # 固定部分
    fixed_start = "GACUA"  # CUGAU的反向互补
    fixed_end = "CUUUC"  # GAAAC的反向互补

    # 非保守序列
    non_conserved_seq = "GAGGCCGAAAGGCC"

    # 生成所有可能的中间序列
    candidates = []

    # 从5'端开始互补
    for length in range(14, -1, -1):
        if length > 0:
            complement_part = reverse_complement(non_conserved_seq[:length])
            # 补齐到14个碱基，用随机碱基填充
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            # 0个碱基互补，全部随机
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    # 从中间开始互补
    mid_point = len(non_conserved_seq) // 2
    for length in range(14, -1, -1):
        if length > 0:
            # 从中间向两边扩展
            start_idx = max(0, mid_point - length // 2)
            end_idx = min(len(non_conserved_seq), start_idx + length)
            complement_part = reverse_complement(non_conserved_seq[start_idx:end_idx])
            # 补齐到14个碱基，用随机碱基填充
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            # 0个碱基互补，全部随机
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    # 从3'端开始互补
    for length in range(14, -1, -1):
        if length > 0:
            complement_part = reverse_complement(non_conserved_seq[-length:])
            # 补齐到14个碱基，用随机碱基填充
            remaining = 14 - length
            random_part = ''.join(random.choice('AUCG') for _ in range(remaining))
            middle = complement_part + random_part
        else:
            # 0个碱基互补，全部随机
            middle = ''.join(random.choice('AUCG') for _ in range(14))

        inhibition_strand = fixed_start + middle + fixed_end
        candidates.append(inhibition_strand)

    return candidates


def calculate_binding_energy(seq1, seq2, model):
    """计算两个序列的结合自由能"""
    # 创建复合物
    complex = Complex([Strand(seq1, name="strand1"), Strand(seq2, name="strand2")])

    # 计算结合自由能
    result = complex_analysis(complexes=[complex], model=model, compute=['pfunc'])

    # 获取自由能（单位：kcal/mol）
    energy = result[complex].free_energy
    return energy


def generate_and_screen_irs(target_rna):
    """生成并筛选IRS序列"""
    # 设置物理模型
    model = Model(material='rna', celsius=37)

    # 生成识别链（目标RNA的反向互补）
    recognition_strand = reverse_complement(target_rna)

    # 生成候选抑制链
    inhibition_candidates = generate_inhibition_strands()

    # HiBiT Switch序列
    hibit_switch = "UCUCCUCUGGCGACCCUGAUGAGGCCGAAAGGCCGAAACGGUAUCGACCGUAGGUUGCCAGAACAGAGGAGAUAAAGAUGGUGAGCGGCUGGCGGCUGUUCAAGAAGAUUAGC"

    # 存储结果
    results = []

    print(f"正在为目标RNA {target_rna} 生成并筛选IRS序列...")

    # 计算每个候选IRS的结合自由能和评分
    for i, inhibition_strand in enumerate(inhibition_candidates):
        # 构建完整IRS序列
        irs_sequence = inhibition_strand + recognition_strand

        # 计算IRS与Target的结合自由能 (ΔG₁)
        dg1 = calculate_binding_energy(irs_sequence, target_rna, model)

        # 计算IRS与HiBiT Switch的结合自由能 (ΔG₂)
        dg2 = calculate_binding_energy(irs_sequence, hibit_switch, model)

        # 计算置换效率得分
        displacement_score = abs(dg1) - abs(dg2)

        # 存储结果
        results.append({
            'irs_sequence': irs_sequence,
            'inhibition_strand': inhibition_strand,
            'recognition_strand': recognition_strand,
            'dg1': dg1,
            'dg2': dg2,
            'displacement_score': displacement_score
        })

        print(f"候选IRS {i + 1}/{len(inhibition_candidates)} 处理完成")

    # 按置换效率得分降序排序
    sorted_results = sorted(results, key=lambda x: x['displacement_score'], reverse=False)

    return sorted_results


def print_results(results, top_n=5):
    """打印筛选结果"""
    print("\n=== IRS序列筛选结果 ===")
    print(f"{'排名':<4} {'IRS序列':<60} {'ΔG(IRS-Target)':<15} {'ΔG(IRS-HiBiT)':<15} {'置换得分':<15}")
    print("-" * 110)

    for i, result in enumerate(results[:top_n], 1):
        print(
            f"{i:<4} {result['irs_sequence']:<60} {result['dg1']:<15.2f} {result['dg2']:<15.2f} {result['displacement_score']:<15.2f}")


# 示例使用
if __name__ == "__main__":
    # 示例目标RNA序列
    target_rna_example = "GAUACUCCUAAUUAUGAUGUGCAGAAACACAUCAA"

    # 生成并筛选IRS序列
    results = generate_and_screen_irs(target_rna_example)

    # 打印结果
    print_results(results, top_n=100)

    # 可选：保存所有结果到文件
    import csv

    with open('irs_screening_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['irs_sequence', 'inhibition_strand', 'recognition_strand',
                      'dg1', 'dg2', 'displacement_score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print("\n所有结果已保存到 irs_screening_results.csv")
