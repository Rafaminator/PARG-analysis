#!/bin/bash
# ============================================================
# Script: find_motifs_exact.sh
# Descrição:
#   Busca motivos específicos em FASTA SEM mismatches,
#   mas todas as buscas subsequentes (1–4, 1–8 e 4–6)
#   são feitas SOMENTE nas reads que contêm o SNP do exon 9.
# ============================================================

set -euo pipefail

echo "=== Iniciando busca por motivos SEM tolerância (100% identidade) ==="

# ------------------------------------------------------------
# Arquivos originais
PARG="reads_PARG.fa"
PARGP1="reads_PARGP1.fa"

# ------------------------------------------------------------
# Pré-limpeza
echo ">> Limpando arquivos FASTA..."
iconv -f utf-8 -t ascii//TRANSLIT -c "$PARG"   | sed 's/[^ACGTNacgtn>]//g' > reads_PARG_clean.fa
iconv -f utf-8 -t ascii//TRANSLIT -c "$PARGP1" | sed 's/[^ACGTNacgtn>]//g' > reads_PARGP1_clean.fa

PARG="reads_PARG_clean.fa"
PARGP1="reads_PARGP1_clean.fa"

# ------------------------------------------------------------
# Motivos
motif_parg="catttccacgacgaaatgctaagatgaaat"
motif_pargp1="catttccacgatgaaatgctcagatgaaat"

motif_14="gtccatctctctcagaaaagaaca"
motif_14_55="gtccatctctctcagataagaagt"
motif_18="acctcgcttggtacttgaag"

# Junções 4–6 (3 variantes)
j46_parg="ACTATTCGGAATGGTGAG"
j46_pseudo="ACTGTTTGGAATGGTGAG"
j46_meyer="ACTATTTGGAATGGTGAG"

# ------------------------------------------------------------
FLAG_EXACT="-i"

search_motif() {
    local motif="$1"
    local fasta="$2"
    local output="$3"
    seqkit grep -s $FLAG_EXACT -p "$motif" "$fasta" > "$output"
}

# ------------------------------------------------------------
echo ">> EXTRAINDO as reads com SNP do exon 9..."

search_motif "$motif_parg"   "$PARG"   "PARG_snp9.fa"
search_motif "$motif_pargp1" "$PARGP1" "PARGP1_snp9.fa"

echo "Reads filtradas:"
echo "PARG_snp9.fa:   $(grep -c '^>' PARG_snp9.fa)"
echo "PARGP1_snp9.fa: $(grep -c '^>' PARGP1_snp9.fa)"

# Agora TODAS as buscas seguintes usam apenas PARG_snp9.fa e PARGP1_snp9.fa
# ------------------------------------------------------------

echo ">> Buscando junção exon 1–4..."
search_motif "$motif_14" "PARG_snp9.fa"   "PARG_snp9_exon14.fa"
search_motif "$motif_14" "PARGP1_snp9.fa" "PARGP1_snp9_exon14.fa"

echo ">> Buscando junção exon 1–4_55..."
search_motif "$motif_14_55" "PARG_snp9.fa"   "PARG_snp9_exon14_55.fa"
search_motif "$motif_14_55" "PARGP1_snp9.fa" "PARGP1_snp9_exon14_55.fa"


echo ">> Buscando junção exon 1–8..."
search_motif "$motif_18" "PARG_snp9.fa"   "PARG_snp9_exon18.fa"
search_motif "$motif_18" "PARGP1_snp9.fa" "PARGP1_snp9_exon18.fa"

echo ">> Buscando junção exon 4–6 (3 variantes)..."

for FASTA in PARG_snp9.fa PARGP1_snp9.fa; do
    prefix=$(echo $FASTA | sed 's/.fa$//')
    search_motif "$j46_parg"   "$FASTA"  "${prefix}_j46_parg.fa"
    search_motif "$j46_pseudo" "$FASTA"  "${prefix}_j46_pseudo.fa"
    search_motif "$j46_meyer"  "$FASTA"  "${prefix}_j46_meyer.fa"
done

# ------------------------------------------------------------
echo
echo "=== Contagem de reads encontrados ==="

echo "--- SNP exon 9 ---"
echo "PARG SNP9:"   $(grep -c '^>' PARG_snp9.fa)
echo "PARGP1 SNP9:" $(grep -c '^>' PARGP1_snp9.fa)

echo "--- Junção 1–4 ---"
echo "PARG 1–4:"     $(grep -c '^>' PARG_snp9_exon14.fa)
echo "PARGP1 1–4:"   $(grep -c '^>' PARGP1_snp9_exon14.fa)

echo "--- Junção 1–4_tipo_55 ---"
echo "PARG 1–4_55:"     $(grep -c '^>' PARG_snp9_exon14_55.fa)
echo "PARGP1 1–4_55:"   $(grep -c '^>' PARGP1_snp9_exon14_55.fa)

echo "--- Junção 1–8 ---"
echo "PARG 1–8:"     $(grep -c '^>' PARG_snp9_exon18.fa)
echo "PARGP1 1–8:"   $(grep -c '^>' PARGP1_snp9_exon18.fa)

echo "--- Junção 4–6 ---"
for gene in PARG_snp9 PARGP1_snp9; do
    echo ">>> $gene"
    echo "parg:"   $(grep -c '^>' ${gene}_j46_parg.fa)
    echo "pseudo:" $(grep -c '^>' ${gene}_j46_pseudo.fa)
    echo "meyer:"  $(grep -c '^>' ${gene}_j46_meyer.fa)
done

echo "=== Busca finalizada ==="

echo
echo "=== Exportando headers das reads classificadas ==="

extract_headers() {
    local input="$1"
    local output="$2"
    if [[ -s "$input" ]]; then
        grep '^>' "$input" | sed 's/^>//' > "$output"
    else
        echo "" > "$output"   # cria arquivo vazio se não houver reads
    fi
    echo "Gerado: $output"
}

# SNP exon 9
extract_headers "PARG_snp9.fa"        "PARG_snp9.headers.txt"
extract_headers "PARGP1_snp9.fa"      "PARGP1_snp9.headers.txt"

# Junção 1–4
extract_headers "PARG_snp9_exon14.fa"     "PARG_exon14.headers.txt"
extract_headers "PARGP1_snp9_exon14.fa"   "PARGP1_exon14.headers.txt"

# Junção 1–8
extract_headers "PARG_snp9_exon18.fa"     "PARG_exon18.headers.txt"
extract_headers "PARGP1_snp9_exon18.fa"   "PARGP1_exon18.headers.txt"

# Junções 4–6 (3 variantes)
for gene in PARG_snp9 PARGP1_snp9; do
    extract_headers "${gene}_j46_parg.fa"   "${gene}_j46_parg.headers.txt"
    extract_headers "${gene}_j46_pseudo.fa" "${gene}_j46_pseudo.headers.txt"
    extract_headers "${gene}_j46_meyer.fa"  "${gene}_j46_meyer.headers.txt"
done

echo "=== Headers exportados no diretório atual ==="

