# Signals & Systems: DCT-Based Image Compression

## Project Purpose

This project implements a **Discrete Cosine Transform (DCT)-based image compression** algorithm similar to JPEG compression. It processes grayscale images by applying DCT, quantization, and Huffman encoding to achieve lossy compression while analyzing signal processing metrics and quality degradation.

## How It Works

The compression pipeline consists of five main stages:

### 1. **Blocking & DCT (Discrete Cosine Transform)**
- Divides input image into 8×8 pixel blocks
- Applies DCT to transform spatial domain data to frequency domain
- Concentrates image energy in low-frequency components

### 2. **Quantization (Lossy Compression)**
- Applies a quantization matrix to DCT coefficients
- Reduces precision of coefficients, especially high-frequency components
- This stage introduces irreversible data loss

### 3. **ZigZag Encoding & RLE (Run Length Encoding)**
- Reorders coefficients in zigzag pattern
- Applies DC differential encoding for DC components
- Encodes consecutive zeros using RLE for better entropy encoding
- Converts data to symbol pairs (Run, Size)

### 4. **Huffman Encoding (Entropy Coding)**
- Encodes symbols using variable-length Huffman codes
- DC and AC components use separate Huffman tables
- Produces the final compressed bit stream

### 5. **Decoder & Quality Analysis**
- Reconstructs image through inverse DCT
- Calculates **MSE** (Mean Squared Error) and **PSNR** (Peak Signal-to-Noise Ratio)
- Quantifies compression ratio and quality loss

## Project Structure

```
Project Code/
├── Project.py                 # Main compression algorithm
├── input_images/              # Folder for images to process
├── output_images/             # Folder for compressed/reconstructed images
├── reports/                   # Individual analysis reports for each image
├── analysis_report.txt        # Legacy report file
└── README.md                  # This file
```

## Features

- **Batch Processing**: Automatically processes all images in `input_images` folder
- **Multiple Formats**: Supports PNG, JPG, JPEG, BMP, GIF, TIFF
- **Detailed Reports**: Generates individual `.txt` reports for each processed image
- **Visualizations**: 
  - 8×8 Block DCT coefficient tables
  - Quantized coefficient matrices
  - RLE (Run Length Encoding) sample output
- **Quality Metrics**: Computes MSE and PSNR for each image
- **Compression Analysis**: Reports original size, compressed size, and compression ratio

## Usage

### Running the Script

```bash
python Project.py
```

### Process Flow

1. Script scans `input_images/` folder for image files
2. For each image:
   - Loads and converts to grayscale
   - Applies level shifting (subtract 128)
   - Performs DCT transformation on 8×8 blocks
   - Quantizes DCT coefficients
   - Applies RLE and Huffman encoding
   - Reconstructs image through inverse DCT
   - Saves processed image to `output_images/`
   - Generates detailed report in `reports/`
   - Displays visualization tables and statistics

## Output Description

### Processed Images (`output_images/`)
- Filename format: `processed_{original_name}.jpg`
- Reconstructed images after compression and decompression

### Report Files (`reports/`)
- Filename format: `{image_name}_report.txt`
- Contains:
  1. Original image information (pixel count, file size)
  2. Quantization impact (active coefficients, data ratio)
  3. RLE encoding results
  4. Huffman encoding statistics (compressed bit count, compression ratio)
  5. Quality metrics (MSE, PSNR in dB)

## Configuration Parameters

Located in `Project.py`:

- **QUANT_VAL**: Quantization matrix value (currently 10)
  - Lower values = higher quality, larger file
  - Higher values = lower quality, smaller file

- **Image Paths**:
  - `INPUT_PATH`: Source images directory
  - `OUTPUT_PATH`: Processed images directory
  - `REPORTS_PATH`: Analysis reports directory

## Key Formulas

**DCT Transform:**
$$T_{u,v} = C_u C_v \sum_{x=0}^{7} \sum_{y=0}^{7} f(x,y) \cos\left(\frac{(2x+1)u\pi}{16}\right) \cos\left(\frac{(2y+1)v\pi}{16}\right)$$

**MSE (Mean Squared Error):**
$$\text{MSE} = \frac{1}{M \times N} \sum_{i=0}^{M-1} \sum_{j=0}^{N-1} [I_{\text{orig}}(i,j) - I_{\text{recon}}(i,j)]^2$$

**PSNR (Peak Signal-to-Noise Ratio):**
$$\text{PSNR} = 10 \log_{10}\left(\frac{255^2}{\text{MSE}}\right) \text{ dB}$$

## Dependencies

```
numpy
Pillow (PIL)
matplotlib
```

Install with:
```bash
pip install numpy pillow matplotlib
```

## Example Output

When processing an image named `baboon.png`:

**Console Output:**
```
[1/1] Processing: baboon.png

1. GİRİŞ (ORİJİNAL) BİLGİLERİ
 - Orijinal Piksel Verisi (Sıkıştırılmamış): 4,194,304 bit
 - Orijinal Dosya Boyutu: 524,288 Bayt
 - Toplam 8x8 Blok Sayısı: 1,024 adet

2. KUANTİZASYON SONRASI BİT DEĞİŞİMİ (Kayıplı Sıkıştırma)
 - Aktif Veri Oranı: %15.23

4. HUFFMAN KODLAMA SONRASI (Entropi Kodlama)
 - GERÇEK ÇIKTI BİT MİKTARI (Sıkıştırılmış): 524,288 bit
 - GERÇEK Sıkıştırma Oranı (Ratio): 8.00:1

5. KALİTE METRİKLERİ
 - MSE (Ortalama Karesel Hata): 42.50
 - PSNR (Tepe Sinyal Gürültü Oranı): 31.85 dB
```

**Files Generated:**
- `output_images/processed_baboon.jpg` - Reconstructed image
- `reports/baboon_report.txt` - Detailed analysis report

## Notes

- This is a **lossy compression** algorithm; quality degradation is expected
- PSNR values typically range from 25-35 dB (lower = more quality loss)
- The implementation uses Turkish (`Türkçe`) comments and variable names
- Huffman tables and quantization matrices are fixed and not adaptive

## Author & Version

- **Project Type**: Signals & Systems Academic Project
- **Date**: February 2026
- **Language**: Python 3.10+
