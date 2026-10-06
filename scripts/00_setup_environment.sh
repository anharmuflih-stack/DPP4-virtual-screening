#!/bin/bash

echo "Memulai instalasi perangkat lunak untuk Molecular Docking (AutoDock Vina)..."

# 1. Mengecek dan menginstal Homebrew (jika belum ada)
if ! command -v brew &> /dev/null
then
    echo "Homebrew tidak ditemukan. Menginstal Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
else
    echo "Homebrew sudah terinstal."
fi

# 2. Menginstal OpenBabel (untuk konversi format kimia dan pembuatan struktur 3D)
echo "Menginstal OpenBabel..."
brew install open-babel

# 3. Menginstal AutoDock Vina
echo "Menginstal AutoDock Vina..."
brew install autodock-vina

# 4. Menginstal Meeko (Pustaka Python resmi dari tim AutoDock untuk menggantikan MGLTools)
# Meeko sangat modern dan jauh lebih mudah digunakan pada Mac untuk menyiapkan file .pdbqt
echo "Menginstal Meeko dan dependensi Python..."
/usr/local/bin/python3 -m pip install meeko scipy numpy

echo "=========================================================="
echo "Instalasi Selesai!"
echo "Anda sekarang memiliki:"
echo "1. openbabel (obabel) : Untuk menghasilkan struktur 3D ligan"
echo "2. vina               : Engine untuk Molecular Docking"
echo "3. meeko              : Tools persiapan protein/ligan (.pdbqt)"
echo "=========================================================="
