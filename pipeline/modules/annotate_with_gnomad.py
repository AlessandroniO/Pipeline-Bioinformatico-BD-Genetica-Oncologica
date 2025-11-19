#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import csv
import subprocess
import shlex

def norm_chrom(chrom: str) -> str:
    """
    Normaliza el nombre del cromosoma para consultar gnomAD.
    - chr23 -> chrX
    - deja todo lo demás igual (chr1, chr2, chr19, etc.)
    """
    if chrom.strip().lower() == "chr23":
        return "chrX"
    return chrom.strip()

def get_gnomad_af(gnomad_vcf, chrom, pos, ref, alt):
    """
    Consulta una coordenada puntual en gnomAD usando bcftools view.
    Devuelve la AF (frecuencia alélica) para el alelo ALT específico, o None si no encuentra.
    """
    # Construimos región tipo chr19:19147134-19147134
    region = f"{chrom}:{pos}-{pos}"

    # bcftools view -r chr19:19147134-19147134 gnomad.vcf.gz
    cmd = f"bcftools view -r {region} {shlex.quote(gnomad_vcf)}"
    # Ejecutamos bcftools inside this Python, capturando stdout
    try:
        res = subprocess.run(
            shlex.split(cmd),
            capture_output=True,
            text=True,
            check=False
        )
    except Exception as e:
        # Si bcftools falla duro
        return None

    if res.returncode != 0:
        # No se pudo leer esa región, o no existe
        return None

    # Parsear líneas VCF devueltas (puede haber varias variantes en misma posición).
    # Nos quedamos solo con la que matchee REF y ALT.
    for line in res.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 8:
            continue
        v_chrom, v_pos, v_id, v_ref, v_alt, v_qual, v_filter, v_info = fields[:8]

        if v_chrom != chrom:
            continue
        if v_pos != str(pos):
            continue
        if v_ref != ref:
            continue

        # Ojo: gnomAD puede tener múltiples ALT separados por coma.
        # Tenemos que identificar cuál ALT coincide y extraer la AF correspondiente.
        alts = v_alt.split(",")
        if alt not in alts:
            continue

        alt_index = alts.index(alt)  # índice 0,1,2,...

        # En gnomAD suele haber AC, AN, AF, etc. AF puede ser coma-separada en el mismo orden que los ALT.
        # Ejemplo de INFO: "AC=5;AF=0.00002;AN=251000;..."
        info_dict = {}
        for kv in v_info.split(";"):
            if "=" in kv:
                k,v = kv.split("=",1)
                info_dict[k] = v
            else:
                info_dict[kv] = True

        if "AF" in info_dict:
            af_values = info_dict["AF"].split(",")
            if alt_index < len(af_values):
                return af_values[alt_index]

        # Si no hay AF en INFO o no coincide índice → no devolvemos nada
    return None

def main(tsv_in, gnomad_vcf_gz, tsv_out):
    """
    Lee somatic_summary.tsv o germline_summary.tsv
    Columnas esperadas (sin header, como las tuyas):
    0 CHROM
    1 POS
    2 REF
    3 ALT
    4 SAMPLE_COUNT
    5 CLINVAR_MATCH
    6 CLINVAR_INFO
    7 SOURCE
    """
    rows_in = []
    with open(tsv_in, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for cols in reader:
            if not cols:
                continue
            # Skip líneas vacías
            if len(cols) < 8:
                # si hay líneas raras, las ignoramos
                continue
            chrom, pos, ref, alt, scount, cmatch, cinfo, source = cols[:8]
            rows_in.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "SAMPLE_COUNT": scount,
                "CLINVAR_MATCH": cmatch,
                "CLINVAR_INFO": cinfo,
                "SOURCE": source
            })

    # Para cada fila, buscamos AF en gnomAD
    annotated = []
    for r in rows_in:
        chrom_norm = norm_chrom(r["CHROM"])
        af = get_gnomad_af(
            gnomad_vcf=gnomad_vcf_gz,
            chrom=chrom_norm,
            pos=r["POS"],
            ref=r["REF"],
            alt=r["ALT"]
        )
        newrow = dict(r)
        newrow["gnomAD_AF"] = af if af is not None else ""
        annotated.append(newrow)

    # Guardar TSV con header
    with open(tsv_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "CHROM","POS","REF","ALT",
            "SAMPLE_COUNT","CLINVAR_MATCH","CLINVAR_INFO","SOURCE",
            "gnomAD_AF"
        ])
        for r in annotated:
            writer.writerow([
                r["CHROM"],
                r["POS"],
                r["REF"],
                r["ALT"],
                r["SAMPLE_COUNT"],
                r["CLINVAR_MATCH"],
                r["CLINVAR_INFO"],
                r["SOURCE"],
                r["gnomAD_AF"]
            ])

    print(f"Listo -> {tsv_out}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python annotate_with_gnomad.py <summary.tsv> <gnomad.vcf.gz> <salida.tsv>")
        sys.exit(1)

    tsv_in = sys.argv[1]
    gnomad_vcf_gz = sys.argv[2]
    tsv_out = sys.argv[3]
    main(tsv_in, gnomad_vcf_gz, tsv_out)
