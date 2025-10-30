input_file = "GRCh38_latest_clinvar.vcf.gz.vcf"
output_file = "GRCh38_subset.vcf"
max_variants = 10000  # cambia este número si querés más o menos

with open(input_file, "r", encoding="utf-8") as fin, open(output_file, "w", encoding="utf-8") as fout:
    count = 0
    for line in fin:
        if line.startswith("##"):
            fout.write(line)
        elif line.startswith("#CHROM"):
            fout.write(line)
        else:
            if count < max_variants:
                fout.write(line)
                count += 1
            else:
                break

print(f"✅ Archivo reducido guardado como {output_file} ({count} variantes)")
