import numpy as np
from PIL import Image
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ====================================================
# BÖLÜM 0: KONFİGÜRASYON VE VERİ YAPILARI
# ====================================================

INPUT_PATH = r"C:\Users\Berkeley\Desktop\Signals&Systems Project\Project Code\input_images"
OUTPUT_PATH = r"C:\Users\Berkeley\Desktop\Signals&Systems Project\Project Code\output_images"
REPORTS_PATH = r"C:\Users\Berkeley\Desktop\Signals&Systems Project\Project Code\reports"
LOG_PATH = r"C:\Users\Berkeley\Desktop\Signals&Systems Project\Project Code\analysis_report.txt"

log_data = {}

def function(IMAGE_PATH):
# --- KUANTİZASYON MATRİSİ ---
    QUANT_VAL = 10
    LUMINANCE_QUANT_MATRIX = np.full((8, 8), QUANT_VAL, dtype=np.float64)

    # --- HUFFMAN KOD TABLOLARI (Luminance) ---
    L_DC_HUFFMAN_TABLE = {
        0: ('0000', 4), 1: ('001', 3), 2: ('010', 3), 3: ('011', 3), 4: ('100', 3),
        5: ('101', 3), 6: ('110', 3), 7: ('1110', 4), 8: ('11110', 5), 9: ('111110', 6),
        10: ('1111110', 7), 11: ('11111110', 8)
    }

    L_AC_HUFFMAN_TABLE = {
        0x00: ('1010', 4), 0x01: ('00', 2), 0x02: ('01', 2), 0x03: ('100', 3), 0x04: ('1011', 4),
        0x05: ('1100', 4), 0x06: ('11011', 5), 0x07: ('11100', 5), 0x08: ('11101', 5), 
        0x09: ('111100', 6), 0x0A: ('1111010', 7), 0xF0: ('11111111001', 11)
    }

    # ====================================================
    # BÖLÜM 1: GİRİŞ VE ANALİZ BAŞLANGICI
    # ====================================================

    try:
        img = Image.open(IMAGE_PATH).convert('L')
        gray_matris = np.array(img).astype(np.float64)
        yukseklik_orijinal, genislik_orijinal = gray_matris.shape
        
        # Rapor verileri 1
        log_data['orijinal_piksel_bit'] = yukseklik_orijinal * genislik_orijinal * 8
        log_data['orijinal_dosya_boyutu'] = os.path.getsize(IMAGE_PATH)
        
        level_shifted_matris = gray_matris - 128
        
    except FileNotFoundError:
        print(f"HATA: Dosya bulunamadı: {IMAGE_PATH}")
        sys.exit(1)

    # ====================================================
    # BÖLÜM 2: BLOKLAMA, DCT VE KUANTİZASYON
    # ====================================================

    def padding_ve_bloklara_ayir(matris):
        y, x = matris.shape
        y_yeni, x_yeni = int(np.ceil(y / 8) * 8), int(np.ceil(x / 8) * 8)
        padded = np.pad(matris, ((0, y_yeni - y), (0, x_yeni - x)), mode='edge')
        
        bloklar = []
        for row in range(0, y_yeni, 8):
            for col in range(0, x_yeni, 8):
                bloklar.append(padded[row:row+8, col:col+8])
        return bloklar, y_yeni, x_yeni

    blok_listesi, yeni_y, yeni_x = padding_ve_bloklara_ayir(level_shifted_matris)
    log_data['toplam_blok_sayisi'] = len(blok_listesi)

    def create_dct_matrix():
        T = np.zeros((8, 8))
        for i in range(8):
            c = np.sqrt(1/8) if i == 0 else np.sqrt(2/8)
            for j in range(8):
                T[i, j] = c * np.cos((2*j+1)*i*np.pi / 16)
        return T

    T_MAT = create_dct_matrix()
    quantized_bloklar = []
    non_zero_coeffs = 0

    for blok in blok_listesi:
        dct = T_MAT.dot(blok).dot(T_MAT.T)
        quant = np.round(dct / LUMINANCE_QUANT_MATRIX).astype(np.int32)
        quantized_bloklar.append(quant)
        non_zero_coeffs += np.count_nonzero(quant)

    # Rapor verileri 2
    log_data['toplam_katsayi_giris'] = len(blok_listesi) * 64
    log_data['aktif_katsayi_sayisi'] = non_zero_coeffs
    log_data['aktif_veri_orani'] = (non_zero_coeffs / log_data['toplam_katsayi_giris']) * 100

    # ====================================================
    # BÖLÜM 3: ZIGZAG VE KODLAMA (RLE/DC)
    # ====================================================

    ZIGZAG = [0, 1, 5, 6, 14, 15, 27, 28, 2, 4, 7, 13, 16, 26, 29, 42, 3, 8, 12, 17, 25, 30, 41, 43, 9, 11, 18, 24, 31, 40, 44, 53, 10, 19, 23, 32, 39, 45, 52, 54, 20, 22, 33, 38, 46, 51, 55, 60, 21, 34, 37, 47, 50, 56, 59, 61, 35, 36, 48, 49, 57, 58, 62, 63]

    def get_cat(v): return 0 if v == 0 else int(np.floor(np.log2(np.abs(v)))) + 1
    def get_sup(v):
        if v == 0: return ""
        b = bin(np.abs(v))[2:]
        return b if v > 0 else "".join('1' if x == '0' else '0' for x in b.zfill(get_cat(v)))

    prev_dc = 0
    huff_stream = []
    symbol_count = 0

    for q in quantized_bloklar:
        zz = q.flatten()[ZIGZAG]
        dc_diff = int(zz[0]) - prev_dc
        prev_dc = int(zz[0])
        
        ac_list, run = [], 0
        for val in zz[1:]:
            if val == 0:
                run += 1
                if run == 16: ac_list.append((0xF0, "")); run = 0
            else:
                ac_list.append(((run << 4) | get_cat(val), get_sup(val))); run = 0
        ac_list.append((0x00, ""))
        
        huff_stream.append({'dc_k': get_cat(dc_diff), 'dc_s': get_sup(dc_diff), 'ac': ac_list})
        symbol_count += 1 + len(ac_list)

    log_data['huffman_giris_sembol'] = symbol_count

    # ====================================================
    # BÖLÜM 4: HUFFMAN VE BİT SAYIMI
    # ====================================================

    final_bits = ""
    for b in huff_stream:
        if b['dc_k'] in L_DC_HUFFMAN_TABLE:
            final_bits += L_DC_HUFFMAN_TABLE[b['dc_k']][0] + b['dc_s']
        for s, sup in b['ac']:
            if s in L_AC_HUFFMAN_TABLE:
                final_bits += L_AC_HUFFMAN_TABLE[s][0] + sup

    log_data['gercek_cikti_bit'] = len(final_bits)
    log_data['sikistirma_orani'] = log_data['orijinal_piksel_bit'] / len(final_bits)

    # ====================================================
    # BÖLÜM 5: DECODER VE KALİTE (MSE/PSNR)
    # ====================================================

    recon = np.zeros((yeni_y, yeni_x))
    for i, q in enumerate(quantized_bloklar):
        r, c = (i // (yeni_x // 8)) * 8, (i % (yeni_x // 8)) * 8
        dq = q * LUMINANCE_QUANT_MATRIX
        recon[r:r+8, c:c+8] = T_MAT.T.dot(dq).dot(T_MAT)

    final_img = np.clip(recon[:yukseklik_orijinal, :genislik_orijinal] + 128, 0, 255).astype(np.uint8)
    output_file = os.path.join(OUTPUT_PATH, f"processed_{os.path.splitext(os.path.basename(IMAGE_PATH))[0]}.jpg")
    Image.fromarray(final_img).save(output_file)

    mse = np.mean((gray_matris - final_img.astype(np.float64))**2)
    psnr = 10 * np.log10(255**2 / mse)

    # ====================================================
    # RAPOR OLUŞTURMA (GÖRSELDEKİ FORMATTA)
    # ====================================================

    def generate_report(d):
        r =  "1. GİRİŞ (ORİJİNAL) BİLGİLERİ\n"
        r += f" - Orijinal Piksel Verisi (Sıkıştırılmamış): {d['orijinal_piksel_bit']:,} bit\n"
        r += f" - Orijinal Dosya Boyutu: {d['orijinal_dosya_boyutu']:,} Bayt\n"
        r += f" - Toplam 8x8 Blok Sayısı: {d['toplam_blok_sayisi']:,} adet\n\n"
        
        r += "2. KUANTİZASYON SONRASI BİT DEĞİŞİMİ (Kayıplı Sıkıştırma)\n"
        r += " Bu aşamada kayıp yaşanır ve verinin büyük kısmı sıfırlanır.\n"
        r += f" - Toplam Katsayı Sayısı (Giriş): {d['toplam_katsayi_giris']:,} adet\n"
        r += f" - **Sıfır Olmayan Katsayı Sayısı (Aktif Veri): {d['aktif_katsayi_sayisi']:,} adet**\n"
        r += f" - Aktif Veri Oranı: %{d['aktif_veri_orani']:.2f}\n\n"
        
        r += "3. RLE VE DC FARK KODLAMA SONRASI\n"
        r += " Sıfır koşuları kısaltılır ve veriler Huffman için sembol çiftlerine dönüştürülür.\n"
        r += f" - Toplam DC/AC Sembol Çifti Sayısı (Huffman Girişi): {d['huffman_giris_sembol']:,} adet\n\n"
        
        r += "4. HUFFMAN KODLAMA SONRASI (Entropi Kodlama)\n"
        r += " Değişken uzunluklu kodlar kullanılarak semboller bit akışına çevrilir.\n"
        r += f" - **GERÇEK ÇIKTI BİT MİKTARI (Sıkıştırılmış): {d['gercek_cikti_bit']:,} bit**\n"
        r += f" - **GERÇEK Sıkıştırma Oranı (Ratio): {d['sikistirma_orani']:.2f}:1**\n\n"
        
        r += "5. KALİTE METRİKLERİ\n"
        r += f" - MSE (Ortalama Karesel Hata): {mse:.2f}\n"
        r += f" - PSNR (Tepe Sinyal Gürültü Oranı): {psnr:.2f} dB\n"
        return r

    final_report = generate_report(log_data)
    print(final_report)
    os.makedirs(REPORTS_PATH, exist_ok=True)
    report_file = os.path.join(REPORTS_PATH, f"{os.path.splitext(os.path.basename(IMAGE_PATH))[0]}_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f: f.write(final_report)

    # =============================================================================================
    # EKLEME: SADECE İLK 4 BLOK İÇİN SAYISAL GÖRSELLEŞTİRME (SİYAH-BEYAZ TABLO)
    # =============================================================================================
    print("\n--- İlk 4 Blok İçin Tablo Oluşturuluyor ---")

    num_blocks_to_show = 4
    selected_blocks = blok_listesi[:num_blocks_to_show] 
    dct_results = []

    for blok in selected_blocks:
        val = T_MAT.dot(blok).dot(T_MAT.T)
        dct_results.append(val)

    final_vis_matrix = np.hstack(dct_results)
    h, w = final_vis_matrix.shape
    plt.figure(figsize=(16, 5))
    ax = plt.gca()
    cmap_white = mcolors.ListedColormap(['white'])
    ax.imshow(np.zeros_like(final_vis_matrix), cmap=cmap_white, interpolation='nearest')

    for r in range(h):
        for c in range(w):
            val = final_vis_matrix[r, c]
            ax.text(c, r, f"{val:.1f}", ha="center", va="center", color="black", fontsize=9, fontweight='normal')

    ax.set_xticks(np.arange(-0.5, w, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, h, 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    for x_grid in range(0, w + 1, 8):
        ax.axvline(x_grid - 0.5, color='black', linewidth=2)

    ax.axhline(-0.5, color='black', linewidth=2)
    ax.axhline(h - 0.5, color='black', linewidth=2)

    ax.set_xticks(np.arange(3.5, w, 8))
    ax.set_xticklabels([f"BLOK {i+1}" for i in range(num_blocks_to_show)], fontsize=12, fontweight='bold')
    ax.xaxis.tick_top()
    ax.set_yticks([])
    ax.tick_params(axis='both', which='both', length=0)
    plt.tight_layout()
    print("İlk 4 bloğun sayısal değerleri gösteriliyor.")
    plt.show()

    # =============================================================================================
    # EKLEME 2: İLK 4 BLOK İÇİN QUANTIZED (KUANTİZE EDİLMİŞ) MATRİS GÖRSELLEŞTİRME
    # =============================================================================================
    print("\n--- İlk 4 Blok İçin QUANTIZED Tablo Oluşturuluyor ---")

    selected_quant_blocks = quantized_bloklar[:4]
    final_quant_vis_matrix = np.hstack(selected_quant_blocks)

    h, w = final_quant_vis_matrix.shape
    plt.figure(figsize=(16, 5)) 
    ax = plt.gca()
    cmap_white = mcolors.ListedColormap(['white'])
    ax.imshow(np.zeros_like(final_quant_vis_matrix), cmap=cmap_white, interpolation='nearest')

    for r in range(h):
        for c in range(w):
            val = final_quant_vis_matrix[r, c]
            ax.text(c, r, f"{int(val)}", ha="center", va="center", color="black", fontsize=10, fontweight='normal')

    ax.set_xticks(np.arange(-0.5, w, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, h, 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    for x_grid in range(0, w + 1, 8):
        ax.axvline(x_grid - 0.5, color='black', linewidth=2)

    ax.axhline(-0.5, color='black', linewidth=2)
    ax.axhline(h - 0.5, color='black', linewidth=2)

    ax.set_title("QUANTIZED Coefficients (Kuantize Edilmiş Veri) - İlk 4 Blok", fontsize=14, pad=20)
    ax.set_xticks(np.arange(3.5, w, 8))
    ax.set_xticklabels([f"BLOK {i+1}" for i in range(4)], fontsize=12, fontweight='bold')
    ax.xaxis.tick_top() 
    ax.set_yticks([])   
    ax.tick_params(axis='both', which='both', length=0)
    plt.tight_layout()
    print("İlk 4 bloğun QUANTIZED değerleri gösteriliyor.")
    plt.show()

    # =============================================================================================
    # EKLEME 3: RLE (RUN LENGTH ENCODING) SONUÇLARINI TEK SATIRDA GÖSTERME
    # =============================================================================================
    print("\n" + "="*80)
    print(" RLE (RUN LENGTH ENCODING) SAMPLE OUTPUT (SINGLE LINE FORMAT)")
    print("Format: (Run, Size) -> Run: Number of skipped zeros, Size: Value category (bits)")
    print(" Note: (0,0) indicates EOB (End of Block), meaning the end of the block.")
    print("="*80)

    # Sadece ilk bloğu (indeks 0) örnek olarak gösteriyoruz.
    ornek_blok_idx = 0
    blok_data = huff_stream[ornek_blok_idx]

    print(f"--- BLOK {ornek_blok_idx + 1} RLE STREAM ---")

    # (Run, Size) çiftlerini bir listede toplayalım
    rle_pairs = []

    for symbol_byte, bit_string in blok_data['ac']:
        # Sembol baytını (byte) çözümle: Üst 4 bit Run, Alt 4 bit Kategori (Size)
        run_len = symbol_byte >> 4  # İlk 4 bit
        size_cat = symbol_byte & 0x0F # Son 4 bit
        
        # Listeye ekle: (Run, Size)
        rle_pairs.append(f"({run_len},{size_cat})")

    # Listeyi virgülle birleştirip tek satırda yazdır
    print(", ".join(rle_pairs))
    print("\n" + "="*80)

    # =============================================================================================
    # EKLEME 4: INPUT IMAGES KLASÖRÜNDEN RESİMLERİ OKUYUP OUTPUT IMAGES KLASÖRÜNE KAYDETME
    # =============================================================================================
    print("\n" + "="*80)
    print(" İNPUT RESİMLERİ İŞLEME VE ÇIKTI KAYDETME")
    print("="*80)

    # Output klasörünü oluştur
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    print(f"\nÇıkış klasörü hazırlandı: {OUTPUT_PATH}")

    # Input klasöründen tüm resimleri al
    if os.path.exists(INPUT_PATH):
        image_files = [f for f in os.listdir(INPUT_PATH) 
                    if f.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif', '.tiff'))]
        
        if image_files:
            print(f"\nBulunan resim sayısı: {len(image_files)}")
            print("-" * 80)
            
            for idx, filename in enumerate(image_files, 1):
                input_file = os.path.join(INPUT_PATH, filename)
                output_file = os.path.join(OUTPUT_PATH, f"processed_{filename.split('.')[0]}.jpg")
                
                try:
                    # Resmi işle
                    print(f"\n[{idx}/{len(image_files)}] İşleniyor: {filename}")
                    
                    img = Image.open(input_file).convert('L')
                    gray_data = np.array(img).astype(np.float64)
                    
                    # Level shift
                    level_shifted = gray_data - 128
                    
                    # Bloklama
                    blok_liste, new_y, new_x = padding_ve_bloklara_ayir(level_shifted)
                    
                    # DCT ve Kuantizasyon
                    quant_liste = []
                    for blok in blok_liste:
                        dct = T_MAT.dot(blok).dot(T_MAT.T)
                        quant = np.round(dct / LUMINANCE_QUANT_MATRIX).astype(np.int32)
                        quant_liste.append(quant)
                    
                    # Inverse DCT ve Reconstruct
                    recon = np.zeros((new_y, new_x))
                    for i, q in enumerate(quant_liste):
                        r, c = (i // (new_x // 8)) * 8, (i % (new_x // 8)) * 8
                        dq = q * LUMINANCE_QUANT_MATRIX
                        recon[r:r+8, c:c+8] = T_MAT.T.dot(dq).dot(T_MAT)
                    
                    # Final image
                    h_orig, w_orig = gray_data.shape
                    final_reconstructed = np.clip(recon[:h_orig, :w_orig] + 128, 0, 255).astype(np.uint8)
                    
                    # Kaydet
                    Image.fromarray(final_reconstructed).save(output_file)
                    print(f"✓ Kaydedildi: {os.path.basename(output_file)}")
                    
                except Exception as e:
                    print(f"✗ HATA {filename}: {str(e)}")
            
            print("\n" + "-" * 80)
            print("✓ TÜM RESİMLER İŞLENDİ VE KAYDEDILDI!")
            print(f"Çıkış klasörü: {OUTPUT_PATH}")
        else:
            print(f"UYARI: {INPUT_PATH} klasöründe resim bulunamadı!")
    else:
        print(f"HATA: {INPUT_PATH} klasörü bulunamadı!")

def main():
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    images = [f for f in os.listdir(INPUT_PATH) if f.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif', '.tiff'))]
    for idx, filename in enumerate(images, 1):
        print(f"[{idx}/{len(images)}] Processing: {filename}")
        function(os.path.join(INPUT_PATH, filename))


if __name__ == "__main__":
    main()