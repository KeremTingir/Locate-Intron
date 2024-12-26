# İntron Bulma Algoritması (Kayan Pencere Yaklaşımı)

Bu proje, verilen bir DNA dizisi ve amino asit dizisine göre, bu dizideki intron aralıklarını bulan bir Python kodudur. Kod, kayan pencere (sliding window) yaklaşımını kullanarak mRNA dizisi üzerindeki intron bölgelerini tespit eder.

## Proje Amacı

Bu proje, biyoinformatik alanında yaygın olarak kullanılan intron/ekzon yapılarını tespit etme sürecini otomatikleştirmeyi amaçlamaktadır. Genetik dizilerdeki kodlanmayan bölgeler olan intronların başlangıç ve bitiş pozisyonlarını bulmak, gen ifadesi ve protein sentezi süreçlerini anlamak için önemlidir.

## Nasıl Çalışır?

1.  **DNA -> mRNA Dönüşümü:** Verilen DNA dizisi, transkripsiyon (yazılma) kurallarına göre mRNA dizisine dönüştürülür.
2.  **Kayan Pencere Yaklaşımı:** mRNA dizisi üzerinde 3'lü nükleotit grupları (kodonlar) alınarak, hedef amino asit dizisi ile karşılaştırılır.
3.  **Eşleşme Kontrolü:**
    *   Eğer bir kodon, hedef amino asit dizisindeki sıradaki amino asidi kodluyorsa, bu bölge bir ekzon olarak kabul edilir ve sonraki kodona geçilir.
    *   Eğer bir kodon, hedef amino asit dizisindeki sıradaki amino asidi kodlamıyorsa (eşleşme yoksa), o bölge intron olarak işaretlenir ve kayan pencere bir nükleotit kaydırılarak devam edilir.
4.  **İntron Birleştirme:** Ardışık gelen intron aralıkları birleştirilir.
5. **Sonuç:** Bulunan ve birleştirilen iki intronun başlangıç ve bitiş indeksleri döndürülür. Eğer şartlar sağlanmıyorsa boş liste döndürülür.

## Kodun Yapısı

Proje, tek bir Python dosyası (`locate_intron.py`) içindeki `locate_intron(gene, amino_acids)` fonksiyonunu içerir. Bu fonksiyon, şu adımları gerçekleştirir:

*   **`transcription_table`:** DNA'dan mRNA'ya dönüşüm kurallarını tutar.
*   **`codon_to_amino`:** mRNA kodonlarını amino asitlere eşleştiren bir tablodur.
*   **`locate_intron(gene, amino_acids)`:** DNA dizisini ve amino asit dizisini alır ve intron aralıklarını bulur.

### `locate_intron(gene, amino_acids)` Fonksiyonunun Detaylı Açıklaması:

Bu fonksiyon, temel olarak aşağıdaki adımları gerçekleştirir:

1.  **DNA'dan mRNA'ya Dönüşüm:**

    *   Verilen `gene` listesindeki DNA nükleotitlerini, `transcription_table` sözlüğünü kullanarak mRNA nükleotitlerine dönüştürür.
    *   DNA'daki "A" nükleotiti mRNA'da "U", "T" nükleotiti "A", "C" nükleotiti "G" ve "G" nükleotiti "C" olarak değiştirilir.
    *   Elde edilen mRNA nükleotitleri birleştirilerek tek bir mRNA dizisi oluşturulur.
    
    ```python
    transcription_table = {
        "A": "U", "T": "A", "C": "G", "G": "C"
    }
    mrna = "".join(transcription_table[nuc] for nuc in gene)
    ```

2.  **Kodon-Amino Asit Eşleştirme Tablosu:**

    *   `codon_to_amino` adlı sözlük, mRNA'daki üçlü nükleotit dizileri (kodonlar) ile karşılık gelen amino asitleri eşleştirir.
    *   Bu tablo, kodonların hangi amino asitleri kodladığını belirlemek için kullanılır.
    
    ```python
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
    ```

3.  **Değişkenlerin Tanımlanması:**

    *   `introns`: Bulunan intron aralıklarını tutacak bir listedir.
    *   `current_intron_start`: Geçerli bir intronun başlangıç indeksini saklar.
    *   `mrna_index`: mRNA dizisinde o an hangi pozisyonda (nükleotit) olduğumuzu takip eder.
    *   `amino_index`: `amino_acids` dizisinde o an hangi amino aside baktığımızı takip eder.

    ```python
    introns = []
    current_intron_start = None
    mrna_index = 0          
    amino_index = 0           
    ```

4.  **Kayan Pencere Döngüsü:**

    *   `while` döngüsü, mRNA dizisi boyunca kayarak ilerler.
    *   Her adımda, mRNA dizisinden 3'lü bir nükleotit grubu (kodon) alınır.
    *   Alınan kodon, `codon_to_amino` tablosu kullanılarak karşılık gelen amino aside çevrilmeye çalışılır.
    *   Eğer çevrilen amino asit, `amino_acids` listesindeki sıradaki amino asitle eşleşirse (yani beklenen amino asit ise):
        *   Bir intron başlangıcı kaydedilmişse (`current_intron_start` doluysa), mevcut intron aralığı `introns` listesine eklenir ve `current_intron_start` sıfırlanır.
        *   Sonraki amino aside geçilir (`amino_index` 1 arttırılır) ve mRNA üzerinde 3 nükleotit ilerlenir (`mrna_index` 3 arttırılır).
    *   Eğer amino asitler eşleşmezse (yani farklı bir amino asit bulunursa):
        *   Eğer bir intron başlangıcı kaydedilmemişse (`current_intron_start` boşsa), mevcut mRNA indeksi `current_intron_start` değişkenine kaydedilir (intron başlamış olur).
        *  mRNA üzerinde 1 nükleotit ilerlenir (`mrna_index` 1 arttırılır), amino asit indisi aynı kalır (aynı amino aside bakmaya devam).

    ```python
    while mrna_index <= len(mrna) - 3 and amino_index < len(amino_acids):
        codon = mrna[mrna_index:mrna_index + 3]
        amino_from_codon = codon_to_amino.get(codon, None)
        if amino_from_codon == amino_acids[amino_index]:
            if current_intron_start is not None:
                introns.append((current_intron_start, mrna_index - 1))
                current_intron_start = None
            amino_index += 1
            mrna_index += 3
        else:
            if current_intron_start is None:
                current_intron_start = mrna_index
            mrna_index += 1
    ```

5.  **Kapanmamış İntron Kontrolü:**

    *   Döngü bittiğinde, eğer hala `current_intron_start` değişkeni doluysa (yani açık kalmış bir intron varsa), kalan tüm mRNA dizisi intron olarak kabul edilir ve aralık `introns` listesine eklenir.

    ```python
    if current_intron_start is not None:
        introns.append((current_intron_start, len(mrna) - 1))
    ```

6.  **İntron Birleştirme:**

    *   Bulunan intronlar bitişik veya üst üste biniyor olabilir. Bu durumları ele almak için aşağıdaki işlemler yapılır:
        *   `merged_introns` adında yeni bir liste oluşturulur.
        *   `introns` listesindeki her bir intron aralığı, `merged_introns` listesine eklenir.
        *   Eğer `merged_introns` listesi boş ise, mevcut intron aralığı doğrudan eklenir.
        *   Eğer liste boş değilse ve son eklenen intronun bitiş indeksi mevcut intronun başlangıç indeksinden daha küçükse, mevcut intron aralığı listeye eklenir.
        *   Eğer son eklenen intron ile mevcut intron bitişik veya üst üste biniyorsa, son eklenen intron birleştirilerek güncellenir.

    ```python
    merged_introns = []
    for start, end in introns:
        if not merged_introns or merged_introns[-1][1] < start - 1:
            merged_introns.append((start, end))
        else:
            merged_introns[-1] = (merged_introns[-1][0], end)
    ```

7. **İntron Sayısı Kontrolü:**

 *  Eğer birleştirme sonrası bulunan intron sayısı 2'den farklı ise boş bir liste döndürülür.

    ```python
    if len(merged_introns) != 2:
        return []
    ```

8.  **Sonuç:**
    *   Birleştirilmiş intron aralıkları döndürülür.

    ```python
    return merged_introns
    ```

## Gerekli Kütüphaneler

Bu proje herhangi bir dış kütüphane gerektirmez. Temel Python fonksiyonlarını kullanır.

## Nasıl Kullanılır?

1.  **Dosyayı İndirin:** `locate_intron.py` dosyasını bilgisayarınıza indirin.
2.  **Çalıştırın:** Python yorumlayıcısı kullanarak dosyayı çalıştırın.
3.  **Veri Girişi:** Örnek kullanım bölümünde olduğu gibi `gene` (DNA dizisi) ve `amino_acids` (amino asit dizisi) listelerini güncelleyin.
4.  **Sonuçları Görüntüleyin:** Program, bulunan intron aralıklarını ekrana basacaktır.

```python
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
    introns = find_introns_sliding(gene, amino_acids)

    # Ekrana bas
    print("İntron aralıkları (kayan pencere yaklaşımı):")
    print(introns)
