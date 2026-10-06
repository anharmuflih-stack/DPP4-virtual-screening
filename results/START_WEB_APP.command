#!/bin/bash

# Dapatkan path root dari project berdasarkan lokasi script ini
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd "$DIR"

echo "========================================"
echo " Memulai Server Web UI Dirof DPP-4..."
echo "========================================"
echo "Mohon tunggu sebentar..."

# Matikan proses flask sebelumnya jika ada yang nyangkut di port 8080
lsof -ti:8080 | xargs kill -9 2>/dev/null

# Coba gunakan python3 dari environment, jika gagal coba path absolut
PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    PYTHON_CMD="/usr/local/bin/python3"
fi

# Jalankan Flask Server di background
$PYTHON_CMD scripts/12_flask_app.py &
FLASK_PID=$!

# Tunggu 2 detik agar server siap
sleep 2

# Cek apakah Flask mati seketika (error)
if ! kill -0 $FLASK_PID 2>/dev/null; then
    echo "========================================"
    echo "❌ ERROR: Server gagal berjalan!"
    echo "Pastikan python3 dan flask sudah terinstal."
    echo "========================================"
    read -p "Tekan Enter untuk keluar..."
    exit 1
fi

# Buka otomatis di browser bawaan Mac
open http://127.0.0.1:8080

echo "Website telah dibuka di browser Anda!"
echo "Biarkan jendela terminal ini terbuka selama Anda menggunakan website."
echo "Untuk mematikan website, tekan CTRL+C lalu tutup jendela ini."

# Tunggu proses background selesai
wait $FLASK_PID
