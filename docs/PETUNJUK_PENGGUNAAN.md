# 🧬 DPP-4 Virtual Screening - Petunjuk Penggunaan

Selamat datang di platform *Computational Drug Discovery* untuk penemuan inhibitor reseptor *Dipeptidyl Peptidase-4 (DPP-4)*. Dokumen ini memuat panduan lengkap tentang cara menggunakan *software* dan mengoperasikan *website* simulasi yang telah dibangun menggunakan arsitektur web performa tinggi (Flask + UI TailwindCSS).

🌐 **Akses Website Publik (Live Hugging Face Spaces):** 
*(Segera hadir setelah deploy)*

---

## 1. Memulai Secara Lokal (Sekali Klik) - STANDAR INDUSTRI UI

Jika Anda menggunakan Mac dan ingin segera menguji antarmuka 3D secara luring (lokal) dengan kecepatan instan:

1. Buka aplikasi **Finder**.
2. Masuk ke *folder* `dpp4_project/results/`.
3. Klik ganda (*double-click*) pada file **`START_WEB_APP.command`**.
4. Sebuah terminal kecil akan muncul, dan peramban web (*browser*) Anda akan langsung membuka alamat `http://localhost:8080`.
5. Nikmati antarmuka UI murni yang sangat luwes dan interaktif.
6. Jika sudah selesai, tutup jendela *browser* dan tekan `CTRL+C` di terminal kecil tersebut lalu tutup.

---

## 2. Publikasi ke Hugging Face Spaces (Gratis, Tanpa Kartu Kredit)

Untuk memamerkan *website* berarsitektur canggih ini ke internet, kita menggunakan **Hugging Face Spaces** yang mendukung teknologi Docker dan menyediakan komputasi bertenaga (RAM 16GB) secara gratis:

**Langkah 1: Siapkan Repositori**
1. Buka [**huggingface.co/spaces**](https://huggingface.co/spaces) dan login/buat akun gratis.
2. Klik tombol **Create new Space** di pojok kanan atas.
3. Isi **Space Name** (misal: `Dirof-DPP4-Lab`).
4. Pada bagian **Select the Space SDK**, pilih **Docker** (lalu pilih *Blank*).
5. Pada bagian *Space Hardware*, pilih yang gratis (2 vCPU, 16GB RAM).
6. Klik **Create Space**.

**Langkah 2: Unggah Kode ke Peladen**
Setelah *Space* Anda jadi, buka terminal di komputer Anda, pastikan Anda berada di direktori `dpp4_project`, lalu jalankan dua baris perintah ini:

```bash
git remote add hf https://huggingface.co/spaces/USERNAME_ANDA/NAMA_SPACE_ANDA
git push hf main
```
*(Ganti tautan di atas dengan URL Git yang diberikan oleh layar Hugging Face Anda. Anda akan diminta memasukkan username & token/password).*

Hugging Face akan secara otomatis membaca `Dockerfile` yang telah kami buatkan, menginstal Linux Vina, dan menyalakan *website* HTML orisinal Anda untuk dapat diakses oleh siapapun di seluruh dunia tanpa batas!

---

## 3. Struktur Data Proyek

- `data/01_raw/`: Berisi pangkalan data *raw* bioaktivitas dari ChEMBL.
- `data/03_pdb/`: Struktur protein DPP-4 (3G0B & 5T4B) bersih yang telah disiapkan.
- `data/04_ligands_pdbqt/`: Himpunan molekul siap tambat dalam format PDBQT.
- `results/`: Seluruh catatan Log *docking* dan luaran konformasi molekul.
- `scripts/`: Berisi *backend* Flask dan skrip otomatisasi laboratorium.

*Selamat melakukan eksplorasi dan penemuan obat baru dengan performa web modern!*
