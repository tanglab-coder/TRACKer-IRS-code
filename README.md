This Script is for Automatic Generation and Screening of IRS (Inhibition-Recognition Strand) Sequences  
This script is used to automatically generate and screen IRS (Inhibition-Recognition Strand) sequences. An IRS is an RNA sequence composed of two parts: the **Inhibition Strand** and the **Recognition Strand**. The goal is to design optimal IRS for different target RNA sequences to achieve precise activation control of ribozymes. The specific algorithm is as follows:  

1️⃣ Design Principles  
1.1 Recognition Strand  
- Consists of the **reverse complementary sequence** of the input target RNA sequence.  
- Follows the base pairing rules:  
  - A ↔ U  
  - C ↔ G  
1.2 Inhibition Strand  
- Fixed length**: 24 nucleotides.  
- Used to pair with the conserved structural region of the HH15 ribozyme and inhibit its activity.  
- Construction principles are as follows:  
- Fully complementary** to the conserved region sequence of the HH15 ribozyme:  
    The conserved sequences include CUGAU and GAAAC.
- Partially complementary** to the non-conserved region (located between the conserved regions):  
    The non-conserved sequence is: GAGGCCGAAAGGCC.  
 - The overall inhibition strand must conform to the following structural template:  
    - GUUUC??????????????AUCAG  
- Among them:  
  - "AUCAG" is the reverse complement of CUGAU;  
  - "GUUUC" is the reverse complement of GAAAC;  
  - For the 14 nucleotide positions (denoted by "??????????????") in the middle, complementary nucleotides are generated in the following ways, and candidate sequences are formed by sorting:  
        1. Complement 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, and 0 nucleotides sequentially starting from the 5’ end of the non-conserved sequence (5’GAGGCCGAAAGGCC3’);  
        2. Complement 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, and 0 nucleotides sequentially starting from the middle of the non-conserved sequence;  
        3. Complement 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, and 0 nucleotides sequentially starting from the 3’ end of the non-conserved sequence.  

2️⃣ IRS Construction  
- Complete IRS sequence = Recognition Strand + Inhibition Strand  
- For each target RNA, multiple candidate inhibition strands are generated, and the corresponding recognition strand is concatenated to construct a set of candidate IRS sequences.  

3️⃣ Binding Free Energy Calculation (Using NUPACK)  
For each IRS sequence, the following calculations are performed:  
1. Binding free energy between IRS and Target (ΔG₁)  
   - Reflects the ability of the recognition strand to be successfully displaced.  
2. Binding free energy between IRS and HiBiT Switch (ΔG₂)  
   - Reflects the inhibitory stability of the inhibition strand on the ribozyme structure.
   
4️⃣ Displacement Efficiency Scoring and Sorting  
Calculate the **Displacement Efficiency Score**, defined as:  
Displacement Score = |ΔG₁| - |ΔG₂|  
- Where:  
  - ΔG₁ = Binding free energy between IRS and Target  
  - ΔG₂ = Binding free energy between IRS and HiBiT Switch  
- A higher score indicates that the IRS is more easily displaced by the target RNA, while maintaining good inhibitory effect and responsiveness.
  
5️⃣ Output Content  
All candidate IRS sequences are sorted in descending order of Displacement Score, with the following information output:  
- IRS nucleotide sequence  
- ΔG (IRS–Target)  
- ΔG (IRS–HiBiT Switch)  
- Displacement Score
  
6️⃣ Reference Sequences of Ribozyme and Switch  
- HH15 Ribozyme Sequence:  
  GGCGACCCUGAUGAGGCCGAAAGGCCGAAACGGU  
- HiBiT Switch Sequence:  
  UCUCCUCUGGCGACCCUGAUGAGGCCGAAAGGCCGAAACGGUAUCGACCGUAGGUUGCCAGAACAGAGGAGAUAAAGAUGGUGAGCGGCUGGCGGCUGUUCAAGAAGAUUAGC  
