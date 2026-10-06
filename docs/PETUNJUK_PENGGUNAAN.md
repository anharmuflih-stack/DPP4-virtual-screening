# 🧬 DPP-4 Virtual Screening - Petunjuk Penggunaan

Selamat datang di platform *Computational Drug Discovery* untuk penemuan inhibitor reseptor *Dipeptidyl Peptidase-4 (DPP-4)*. Dokumen ini memuat panduan lengkap tentang cara menggunakan *software* dan mengoperasikan *website* simulasi yang telah dibangun.

🌐 **Akses Website Publik (Live):** 
[**https://dpp4virtualscreener.streamlit.app**](https://dpp4virtualscreener.streamlit.app)

---

## 1. Memulai Secara Lokal (Sekali Klik)

Jika Anda menggunakan Mac dan ingin segera menguji antarmuka 3D secara luring (lokal), kami telah menyediakan *shortcut* khusus yang tidak memerlukan keahlian pemrograman:

1. Buka aplikasi **Finder**.
2. Masuk ke *folder* `dpp4_project/results/`.
3. Klik ganda (*double-click*) pada file **`START_WEB_APP.command`**.
4. Sebuah terminal kecil akan muncul, dan peramban web (*browser*) Anda akan langsung membuka alamat `http://localhost:8080`.
5. Anda bebas melakukan simulasi Docking (memasukkan kode SMILES obat). 
6. Jika sudah selesai, tutup jendela *browser* dan tekan `CTRL+C` di terminal kecil tersebut lalu tutup.

---

## 2. Menjalankan Versi Streamlit (Cloud Ready)

Selain versi Flask (*One-Click*), proyek ini juga dilengkapi dengan antarmuka Streamlit murni yang didesain agar tahan banting di lingkungan awan (*cloud*) seperti *Streamlit Community Cloud*.

Jika Anda ingin menjalankan versi ini di komputer lokal:
1. Buka Terminal Mac Anda.
2. Masuk ke ruang kerja (*workspace*) proyek: 
   ```bash
   cd /Users/aanmuf/drug_discovery/dpp4_project
   ```
3. Nyalakan mesin Streamlit:
   ```bash
   streamlit run scripts/13_streamlit_app.py
   ```
4. Jendela peramban akan otomatis terbuka di port `8501`.

---

## 3. Publikasi ke Streamlit Cloud (Online)

Untuk memamerkan karya Anda ke seluruh dunia, lakukan publikasi ke peladen awan:
1. Pastikan seluruh perubahan kode sudah diunggah (*push*) ke akun GitHub Anda.
2. Kunjungi [**share.streamlit.io**](https://share.streamlit.io/)
3. Klik tombol **New App**.
4. Pilih repositori `dpp4_project` milik Anda.
5. Pada kolom **Main file path**, masukkan tepat teks berikut:
   `scripts/13_streamlit_app.py`
6. Klik **Deploy!**

*Server* Linux pada Streamlit akan mendeteksi kebutuhan proyek, dan mengunduh aplikasi `vina_linux` secara otomatis. Setelah selesai, *website* siap diakses publik!

---

## 4. Struktur Data

- `data/01_raw/`: Berisi pangkalan data *raw* bioaktivitas dari ChEMBL.
- `data/03_pdb/`: Struktur protein DPP-4 (3G0B & 5T4B) bersih yang telah disiapkan.
- `data/04_ligands_pdbqt/`: Himpunan molekul siap tambat dalam format PDBQT.
- `results/`: Seluruh catatan Log *docking* dan luaran konformasi molekul.

*Selamat melakukan eksplorasi dan penemuan obat baru!*
