#!/bin/bash
# ============================================================
# Script: quantify_PARG_isoforms_exact.sh
# Descrição:
#   Quantifica a proporção de reads oriundos do gene PARG
#   que contêm junções exon-exon específicas de cada isoforma.
#   Usa apenas busca exata (sem mismatches, 100% identidade).
# ============================================================

set -euo pipefail

echo "=== Iniciando quantificação das isoformas do gene PARG (modo exato) ==="

# ------------------------------------------------------------
# Arquivo de entrada (gene PARG)
PARG="reads_PARG.fa"

# ------------------------------------------------------------
# Limpeza de caracteres inválidos
echo ">> Limpando arquivo FASTA para evitar erros do seqkit..."
iconv -f utf-8 -t ascii//TRANSLIT -c "$PARG" | sed 's/[^ACGTNacgtn>@]//g' > reads_PARG_clean.fa
PARG="reads_PARG_clean.fa"

# ------------------------------------------------------------
# Motivos de junção exon-exon (confirmados)
motif_PARG111="gccacctcgcttgttttcaaaca"        # junção e1–e2
motif_PARG102a="cattgaggcagttttcaaaca"         # junção e1–e2 (variante 1)
motif_PARG102b="atctctctcaggttttcaaac"         # junção e1–e2 (variante 2)
motif_PARG99="gccacctcgcttgtttggatagtaaag"     # junção e1–e3
motif_PARG53="acctcgcttggtacttgaag"            # junção e1–e8 nova

# ------------------------------------------------------------
# Flag de busca (exata, sem tolerância)
FLAG_EXACT="-i"

# ------------------------------------------------------------
# Função de busca exata
search_motif() {
    local motif="$1"
    local fasta="$2"
    local output="$3"

    echo "🔍 Buscando motivo (exato): $motif"
    seqkit grep -s $FLAG_EXACT -p "$motif" "$fasta" > "$output"
    echo "resposta obtida no modo: exato (-i)"
}

# ------------------------------------------------------------
# Busca por cada isoforma
echo ">> Buscando junções específicas das isoformas..."
search_motif "$motif_PARG111" "$PARG" "PARG111.fa"
search_motif "$motif_PARG102a" "$PARG" "PARG102a.fa"
search_motif "$motif_PARG102b" "$PARG" "PARG102b.fa"
search_motif "$motif_PARG99"   "$PARG" "PARG99.fa"
search_motif "$motif_PARG53"   "$PARG" "PARG53.fa"

# ------------------------------------------------------------
# Contagem de reads
echo
echo ">> Contando reads por isoforma..."
total=$(grep -c "^>" "$PARG")
n111=$(grep -c "^>" PARG111.fa)
n102a=$(grep -c "^>" PARG102a.fa)
n102b=$(grep -c "^>" PARG102b.fa)
n99=$(grep -c "^>" PARG99.fa)
n53=$(grep -c "^>" PARG53.fa)
n102=$((n102a + n102b))

# ------------------------------------------------------------
# Cálculo de percentuais
p111=$(awk -v n="$n111" -v t="$total" 'BEGIN{printf "%.2f", (n/t)*100}')
p102=$(awk -v n="$n102" -v t="$total" 'BEGIN{printf "%.2f", (n/t)*100}')
p99=$(awk -v n="$n99"  -v t="$total" 'BEGIN{printf "%.2f", (n/t)*100}')
p53=$(awk -v n="$n53"  -v t="$total" 'BEGIN{printf "%.2f", (n/t)*100}')

# ------------------------------------------------------------
# Relatório final
echo
echo "=== Quantificação das isoformas do gene PARG (modo exato) ==="
echo "Total de reads analisados (baixa confiança): $total"
echo
echo "PARG111 (e1–e2):             $n111 reads  (${p111}%)"
echo "PARG102 (e1–e2; soma 102a+b): $n102 reads  (${p102}%)"
echo "PARG99  (e1–e3):              $n99 reads  (${p99}%)"
echo "PARG53  (e1–e8 nova):         $n53 reads  (${p53}%)"
echo
echo "=== Análise concluída com sucesso ==="

