# progressive-msa

[English](README.md) | Türkçe

[![CI](https://github.com/aygunhasanbasri8-ui/progressive-msa/actions/workflows/ci.yml/badge.svg)](https://github.com/aygunhasanbasri8-ui/progressive-msa/actions/workflows/ci.yml)

Clustal Omega'nın genel stratejisini izleyen, bağımlılığı olmayan küçük bir Python kütüphanesi:
çoklu dizi hizalama (progressive multiple sequence alignment, MSA).

1. **İkili uzaklıklar:** global Needleman-Wunsch hizalaması ile.
2. **Rehber ağaç:** UPGMA kümelemesi ile.
3. **Aşamalı hizalama:** rehber ağaç boyunca diziler/profiller hizalanır.

Klasik hattın küçük ve okunaklı bir uygulaması; öğrenmek ve küçük girdiler için tasarlandı.
Gerçek biyolojik veride kullanmadan önce [Sınırlar](#sınırlar) bölümüne bakın.

## Özellikler

- Saf Python, çalışma zamanı bağımlılığı yok, Python 3.9 ve üzeri.
- Deterministik sonuç: hizalamadaki ve UPGMA çift seçimindeki eşitlikler her zaman aynı şekilde çözülür.
- Her adım tek başına da çağrılabilir (`calculate_distance_matrix`, `build_guide_tree`, `progressive_alignment`).
- İsteğe bağlı ilerleme günlüğü.

## Kurulum

```bash
pip install -e .
```

## Kullanım

```python
from progressive_msa import ProgressiveMSA

sequences = {
    "Seq1": "ACTG",
    "Seq2": "ACTC",
    "Seq3": "ACCGT",
}

msa = ProgressiveMSA(sequences)
aligned = msa.run()

for name, seq in aligned.items():
    print(name, seq)
```

Çıktı (`verbose=False` verilmezse önce ilerleme günlüğü basılır):

```
Seq3 ACCGT
Seq1 A-CTG
Seq2 A-CTC
```

İlerleme günlüğünü kapatmak için `ProgressiveMSA(...)` çağrısına `verbose=False` verin.

Dönen sözlük, dizileri girdi sırasında değil rehber ağaçta dolaşılma sırasıyla verir.

### API

| Metot | Açıklama |
|---|---|
| `ProgressiveMSA(sequences, verbose=True)` | `sequences` ad → dizi eşlemesidir |
| `run()` | Üç aşamayı çalıştırır, `{ad: hizalanmış_dizi}` döner |
| `calculate_distance_matrix()` | İkili düzenleme uzaklığı matrisi (simetrik, köşegen sıfır) |
| `needleman_wunsch_distance(seq1, seq2)` | İki dizinin global hizalama düzenleme uzaklığı |
| `build_guide_tree()` | İç içe tuple olarak UPGMA rehber ağacı, ör. `('Seq3', ('Seq1', 'Seq2'))`; `calculate_distance_matrix()`'ten sonra çağrılır |
| `progressive_alignment()` | Rehber ağaç boyunca hizalar; `build_guide_tree()`'den sonra çağrılır |

## Nasıl çalışıyor

1. **Uzaklık matrisi.** Her dizi çifti birim maliyetlerle (eşleşmeme 1, boşluk 1) Needleman-Wunsch ile
   hizalanır; çıkan düzenleme uzaklığı simetrik bir matrise yazılır.
2. **Rehber ağaç.** UPGMA, en yakın iki kümeyi birleştirir ve diğer kümelere olan uzaklıklarının
   ortalamasını alır; bunu tek küme kalana dek tekrarlar.
3. **Aşamalı hizalama.** Ağaç aşağıdan yukarı dolaşılır. Her iç düğümde iki alt profil, konsensus
   dizileri (her sütunun en sık kalıntısı) üzerinden hizalanır ve çıkan boşluk deseni iki profilin
   tüm dizilerine uygulanır.

## Sınırlar

- Skorlama basit: eşleşmeme ve boşluk 1 maliyetli; yer değiştirme matrisi (BLOSUM/PAM) ve affine boşluk cezası yok.
- Profiller sütun sütun değil konsensus dizisi üzerinden hizalanıyor; sonuç gerçek bir profil-profil hizalaması değil.
- Uzaklık matrisi, uzunluğu L olan n dizi için O(n²·L²) sürer; büyük girdiler için uygun değil.
- UPGMA sabit bir evrim hızı varsayar.

## Geliştirme

```bash
pip install -e .
pip install pytest
pytest
```

GitHub Actions her push ve pull request'te testleri Python 3.9, 3.11 ve 3.12 üzerinde çalıştırır
(`.github/workflows/ci.yml`).

## Lisans

MIT. Ayrıntılar için [LICENSE](LICENSE).
