# convert_to_genomic_hgvs.py
import sys, time, re
import pandas as pd

import hgvs.parser as P
import hgvs.dataproviders.uta as uta
from hgvs.variantmapper import VariantMapper


def extract_hgvs_c(raw: str):
    """
    Toma cosas tipo:
      'NM_177438.3(DICER1):c.5540C>T (p.Ala1847Val)'
    y devuelve:
      'NM_177438.3:c.5540C>T'

    También soporta si ya viene como 'NM_177438.3:c.5540C>T'
    """
    if not isinstance(raw, str):
        return None

    # caso con gen entre paréntesis
    m = re.search(r'(N[MR]_\d+\.\d+)\([A-Za-z0-9_-]+\):((?:c|n)\.[^ )]+)', raw)
    if m:
        return f"{m.group(1)}:{m.group(2)}"

    # caso ya limpio NM_xxx:c.xxx
    m2 = re.search(r'(N[MR]_\d+\.\d+):(c\.[^ )]+)', raw)
    if m2:
        return f"{m2.group(1)}:{m2.group(2)}"

    return None


def _nc_version(ac: str) -> int:
    """
    'NC_000017.11' -> 11
    Lo usamos para preferir la versión más alta (que suele corresponder a GRCh38).
    """
    try:
        return int(ac.split(".")[-1])
    except Exception:
        return -1


def pick_grch38_alt_ac(hdp, tx_ac: str) -> str:
    """
    Dado un transcripto (ej 'NM_177438.3'), preguntamos a UTA las opciones
    de mapeo a genoma y elegimos el mejor cromosoma de referencia (alt_ac).

    NOTA: get_tx_mapping_options devuelve DictRow (tipo dict-like), NO objetos.
    """
    rows = list(hdp.get_tx_mapping_options(tx_ac))
    if not rows:
        raise RuntimeError(f"Sin opciones de mapeo para {tx_ac}")

    # quedarnos preferentemente con cromosomas 'NC_...'
    nc_rows = [r for r in rows if str(r.get("alt_ac", "")).startswith("NC_")]
    pool = nc_rows if nc_rows else rows

    # elegimos el de versión más alta
    pool_sorted = sorted(
        pool,
        key=lambda r: _nc_version(str(r.get("alt_ac", ""))),
        reverse=True
    )

    return pool_sorted[0]["alt_ac"]


def convert_row(hp, vm, hdp, raw_variant: str):
    """
    Dada la celda 'variante' original:
    1. extrae NM_xxx:c.xxx
    2. parsea esa c. con hgvs
    3. obtiene el cromosoma (alt_ac) desde UTA
    4. proyecta c. -> g. usando VariantMapper

    Devuelve (hgvs_c, hgvs_g) como strings.
    SIN normalización extra.
    """
    hgvs_c = extract_hgvs_c(raw_variant)
    if not hgvs_c:
        raise ValueError("No pude extraer NM_:c. a partir de 'variante'")

    # parseo HGVS c.
    var_c = hp.parse_hgvs_variant(hgvs_c)

    # elegimos genoma (cromosoma) para mapear
    alt_ac = pick_grch38_alt_ac(hdp, var_c.ac)

    # c. -> g.
    var_g = vm.c_to_g(var_c, alt_ac)

    # devolvemos la string directa
    return hgvs_c, str(var_g)


def main(input_csv, output_csv, sleep_s=0.0):
    df = pd.read_csv(input_csv)

    if "variante" not in df.columns:
        print("❌ El CSV debe tener una columna llamada 'variante'")
        sys.exit(1)

    # inicializamos una sola vez
    hp = P.Parser()
    hdp = uta.connect()
    vm  = VariantMapper(hdp)

    uniques = df["variante"].drop_duplicates().tolist()
    map_hgvs_c = {}
    map_hgvs_g = {}
    errors = []

    for raw in uniques:
        try:
            hgvs_c, hgvs_g = convert_row(hp, vm, hdp, raw)
            map_hgvs_c[raw] = hgvs_c
            map_hgvs_g[raw] = hgvs_g
        except Exception as e:
            # si falla (por CNVs raras, problemas de normalización, etc.)
            map_hgvs_c[raw] = extract_hgvs_c(raw) or ""
            map_hgvs_g[raw] = None
            errors.append(f"{raw}\t{e}")

        if sleep_s:
            time.sleep(sleep_s)

    # agregamos columnas nuevas
    df["HGVS_c_extracted"] = df["variante"].map(map_hgvs_c)
    df["HGVS_G_GRCh38"]    = df["variante"].map(map_hgvs_g)

    # guardamos
    df.to_csv(output_csv, index=False)
    print(f"✅ Guardado: {output_csv}")

    if errors:
        with open("error_log.txt", "a", encoding="utf-8") as log:
            for line in errors:
                log.write(line + "\n")
        print("⚠️ Revisar error_log.txt (hubo variantes que no se pudieron proyectar a g.)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python convert_to_genomic_hgvs.py <input.csv> <output.csv>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
