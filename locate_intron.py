def locate_intron(gene, amino_acids):
    
    # 1) DNA -> mRNA dönüştürme tablosu
    transcription_table = {
        "A": "U", "T": "A", "C": "G", "G": "C"
    }

    # DNA'yı mRNA'ya çevir
    mrna = "".join(transcription_table[nuc] for nuc in gene)

    # 2) Kodon -> aminoasit dönüşüm tablosu
    codon_to_amino = {
        "UUU": "PHE", "UUC": "PHE", "UUA": "LEU", "UUG": "LEU",
        "CUU": "LEU", "CUC": "LEU", "CUA": "LEU", "CUG": "LEU",
        "AUU": "ILE", "AUC": "ILE", "AUA": "ILE", "AUG": "MET",
        "GUU": "VAL", "GUC": "VAL", "GUA": "VAL", "GUG": "VAL",
        "UCU": "SER", "UCC": "SER", "UCA": "SER", "UCG": "SER",
        "CCU": "PRO", "CCC": "PRO", "CCA": "PRO", "CCG": "PRO",
        "ACU": "THR", "ACC": "THR", "ACA": "THR", "ACG": "THR",
        "GCU": "ALA", "GCC": "ALA", "GCA": "ALA", "GCG": "ALA",
        "UAU": "TYR", "UAC": "TYR", "UAA": "stop", "UAG": "stop",
        "CAU": "HIS", "CAC": "HIS", "CAA": "GLN", "CAG": "GLN",
        "AAU": "ASN", "AAC": "ASN", "AAA": "LYS", "AAG": "LYS",
        "GAU": "ASP", "GAC": "ASP", "GAA": "GLU", "GAG": "GLU",
        "UGU": "CYS", "UGC": "CYS", "UGA": "stop", "UGG": "TRP",
        "CGU": "ARG", "CGC": "ARG", "CGA": "ARG", "CGG": "ARG",
        "AGU": "SER", "AGC": "SER", "AGA": "ARG", "AGG": "ARG",
        "GGU": "GLY", "GGC": "GLY", "GGA": "GLY", "GGG": "GLY"
    }

    # 3) Değişkenler
    introns = []
    current_intron_start = None
    mrna_index = 0            # mRNA üzerindeki mevcut index
    amino_index = 0           # amino_asits dizisinde mevcut index

    # 4) Kayan pencere döngüsü
    while mrna_index <= len(mrna) - 3 and amino_index < len(amino_acids):
        codon = mrna[mrna_index:mrna_index + 3]

        # Bu 3'lünün geçerli amino asit karşılığı varsa al
        amino_from_codon = codon_to_amino.get(codon, None)

        if amino_from_codon == amino_acids[amino_index]:
            # Eşleşme bulundu!
            # Eğer intron açıksa kapat
            if current_intron_start is not None:
                introns.append((current_intron_start, mrna_index - 1))
                current_intron_start = None

            # Sonraki amino aside geç, kodonu de 3 ilerlet
            amino_index += 1
            mrna_index += 3

        else:
            # Eşleşme yok => intron
            if current_intron_start is None:
                current_intron_start = mrna_index

            # mRNA'yı 1 baz kaydır (frameshift)
            mrna_index += 1

    # 5) Döngü bitti. Hâlâ intron içindeysek (kapanmadıysa) ekle
    if current_intron_start is not None:
        # Kalan tüm bazları intron sayıyoruz
        introns.append((current_intron_start, len(mrna) - 1))

    # 6) İntron aralıklarını birleştirelim (bitişik veya üst üste binen aralıklar)
    merged_introns = []
    for start, end in introns:
        if not merged_introns or merged_introns[-1][1] < start - 1:
            merged_introns.append((start, end))
        else:
            # Son eklenen intron ile birleşiyor
            merged_introns[-1] = (merged_introns[-1][0], end)

    return merged_introns


#-------------------------------------------------------------------
# Örnek kullanım
if __name__ == "__main__":
    # Örnek DNA
    gene = [
        "T","C","T","G","C","A","G","C","A","G","A","G","G","G","G","C","C",
        "G","T","C","G","G","C","A","G","A","A","G","G","A","G","G","G","C",
        "T","C","G","G","G","C","A","G","G","C","T","C","T","G","C","G","A",
        "C","T","C","G","T","A","G","G","C","A","C","C","A","G","G","C","G",
        "T","G","A","G","A","C","C","T","G","T","A","G","C","C","C","C","C",
        "G","A","T","C","A","C","C","A","T","G","T","A","C","A","G","C","T",
        "T","C","A","T","G","G","G","T","G","G","T","G","G","C","C","T","G",
        "T","T","C","T","G","T","G","C","C","T","G","G","G","T","G","G","G",
        "G","A","C","C","A","T","C","C","T","C","C","T","G","G","T","G","G",
        "T","G","G","C","C","A","T","G","G","C","A","A","C","A","G","A","C",
        "G","G","G","G","C","C","A","A","G","G","A","C","A","C","C","T","G",
        "T","A","T","T","C","C","A","G","A","T","G","G","A","G","A","A","C",
        "T","C","T","G","C","G","G","C","T","C","A","A","A","G","A","G","G",
        "G","A","A","A","G","G","G","A","G","C","A","A","C","C","C","A","A",
        "G","G","T","C","A","C","T","C","A","G","C","G","G","A","G","G","C",
        "T","G","A","C","T","C","C","T","G","G","T","C","C","T","A","G","G",
        "C","T","G","G","A","A","G","G","A","G","G","A","A","G","A","A","T",
        "A","G","G","G","C","C","C","A","T","G","G","G","A","G","G","G","A",
        "G","C","T","G","A","G","A","A","G","A","C","T"
    ]

    # Örnek Amino Asit Dizisi
    amino_acids = [
        "ARG","ARG","ARG","LEU","PRO","GLY","SER","ARG","LEU","PRO","PRO",
        "GLU","PRO","VAL","ARG","ASP","ALA","GLU","LEU","VAL","VAL","HIS",
        "VAL","GLU","VAL","PRO","THR","THR","GLY","GLN","ASP","THR","ASP",
        "PRO","PRO","LEU","VAL","GLY","GLY","PRO","PRO","PRO","VAL","PRO",
        "LEU","SER","PRO","PRO","THR","GLU","ASP","GLN","ASP","PRO","THR",
        "PHE","LEU","LEU","LEU","ILE","PRO","GLY","THR","LEU","PRO","ARG",
        "LEU","PHE","stop"
    ]

    # Fonksiyonu çağır
    introns = locate_intron(gene, amino_acids)

    # Ekrana bas
    print("İntron aralıkları (kayan pencere yaklaşımı):")
    print(introns)
