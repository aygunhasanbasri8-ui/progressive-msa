from basriMSA import basriMSA

def main():
    diziler = {
        'Dizi1': 'ACTG',
        'Dizi2': 'ACTC',
        'Dizi3': 'ACCGT'
    }

    msa = basriMSA(diziler)
    sonuclar = msa.run()

    print("\nSonuçlar:")
    for isim, dizi in sonuclar.items():
        print(f"{isim}: {dizi}")

if __name__ == "__main__":
    main()