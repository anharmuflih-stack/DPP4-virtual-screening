# 🧬 DPP-4 Virtual Screening - Petunjuk Penggunaan

Selamat datang di platform *Computational Drug Discovery* untuk penemuan inhibitor reseptor *Dipeptidyl Peptidase-4 (DPP-4)*. Dokumen ini memuat panduan lengkap tentang cara menggunakan *software* dan mengoperasikan *website* simulasi yang telah dibangun menggunakan mesin **Streamlit**.

🌐 **Akses Website Publik (Streamlit Cloud):** 
[**https://dpp4virtualscreener.streamlit.app**](https://dpp4virtualscreener.streamlit.app)

💻 **Kode Sumber / Repositori GitHub:** 
[**https://github.com/anharmuflih-stack/DPP4-virtual-screening**](https://github.com/anharmuflih-stack/DPP4-virtual-screening)

---

## 1. Menjalankan Aplikasi Secara Lokal (Offline)

Proyek ini menyediakan dua mode antarmuka untuk dijalankan secara lokal di komputer (Mac) Anda: **Mode Flask** (Desain orisinal) dan **Mode Streamlit**.
> **Catatan Penting:** Kedua mode aplikasi ini telah terintegrasi secara *real-time* dengan mesin pelacak **PubChem PUG REST API** (Amerika Serikat) dan **EBI ChEMBL API** (Eropa). Pastikan komputer Anda terhubung dengan internet agar sistem dapat menarik data Nama Obat, PubChem CID, dan ChEMBL ID secara otomatis berdasarkan struktur molekul (*SMILES / InChIKey*) yang Anda masukkan.

### A. Mode Localhost Flask (Sekali Klik)
Pilihan terbaik untuk merasakan antarmuka *Tailwind CSS* yang 100% mulus dan responsif.
1. Buka aplikasi **Finder**.
2. Masuk ke *folder* `dpp4_project/results/`.
3. Klik ganda (*double-click*) pada file **`START_WEB_APP.command`**.
4. Sebuah terminal kecil akan muncul, dan peramban web (*browser*) Anda akan langsung membuka alamat `http://localhost:8080`.
5. Untuk mematikan server, tutup jendela *browser* dan tekan `CTRL+C` di terminal kecil tersebut.

### B. Mode Streamlit Lokal
1. Buka aplikasi **Terminal**.
2. Masuk ke ruang kerja proyek Anda:
   ```bash
   cd /Users/aanmuf/drug_discovery/dpp4_project
   ```
3. Nyalakan mesin Streamlit dengan perintah:
   ```bash
   streamlit run scripts/13_streamlit_app.py
   ```
4. Jendela peramban web (*browser*) Anda akan otomatis terbuka dan menampilkan antarmuka aplikasi di port `8501`.
5. Untuk mematikan server lokal, tekan `CTRL+C` di Terminal.

---

## 2. Publikasi ke Streamlit Community Cloud (Online)

Untuk menjaga *website* Anda tetap mengudara di internet secara permanen dan gratis, ikuti langkah publikasi ke **Streamlit Cloud** berikut ini:

1. Pastikan semua *update* atau perubahan kode terbaru sudah diunggah ke GitHub Anda (gunakan perintah `git push` di terminal).
2. Kunjungi [**share.streamlit.io**](https://share.streamlit.io/) dan login menggunakan akun GitHub Anda.
3. Klik tombol **New App**.
4. Pilih repositori `dpp4_project` milik Anda dari daftar yang tersedia.
5. Pada kolom **Main file path**, masukkan dengan teks berikut:
   `scripts/13_streamlit_app.py`
6. Klik **Deploy!**

Server Streamlit Cloud (berbasis Linux) secara otomatis akan mengunduh aplikasi `vina_linux` melalui *script* yang telah disematkan, menyiapkan RDKit, dan memublikasikan 3D *viewer*-nya untuk bisa diakses publik.

---

## 3. Struktur Direktori Data

- `data/01_raw/`: Berisi pangkalan data *raw* bioaktivitas dari ChEMBL.
- `data/03_pdb/`: Struktur protein DPP-4 (3G0B & 5T4B) bersih yang telah disiapkan.
- `data/04_ligands_pdbqt/`: Himpunan molekul siap tambat dalam format PDBQT.
- `results/`: Seluruh catatan Log *docking* dan luaran konformasi molekul (jika dijalankan dengan metode *batch*).
- `scripts/`: Berisi seluruh mesin aplikasi web (`13_streamlit_app.py`) dan alat-alat otomatisasi lainnya.

*Selamat melakukan eksperimen, semoga penemuan obat Anda berhasil!*
