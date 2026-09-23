# -*- coding: utf-8 -*-
"""
Exporta a area de cada unidade amostral com precisao total
==========================================================

POR QUE ISTO EXISTE

A coluna `AOH_unit_area_ha` de `Data_Raw_FINAL.csv` esta gravada como numero
inteiro. Das 67 unidades, 34 sao multiplos exatos de 10, 28 sao multiplos exatos
de 100 e cinco sao multiplos exatos de 1.000 (27000, 33800, 36000, 51000,
150000). Uma area calculada pelo ArcGIS nao produz esse padrao: ela sai como um
numero de ponto flutuante com muitas casas, do tipo 46748.32714. Os valores
redondos indicam que a area foi arredondada, ou digitada, em algum ponto entre o
ArcGIS e a planilha.

Por que importa neste artigo em particular. O argumento central da Parte V e uma
identidade algebrica entre a resposta e um preditor:

    log(OHE) = log(NP) + log(100) - log(PD)

A tabela S28e reporta que essa identidade fecha com discrepancia maxima de 0,029
em escala logaritmica, ou seja 2,9 por cento. Parte dessa discrepancia e
inevitavel: as metricas de paisagem sao contadas sobre pixels de 30 m e a
resposta e a area do poligono, de modo que as duas nao se referem exatamente a
mesma superficie. A outra parte e o arredondamento da resposta, e essa parte e
evitavel. Reportar a area com precisao total elimina uma das duas fontes, deixa
a verificacao da identidade mais limpa e remove uma pergunta obvia de revisor.

Efeito esperado sobre os resultados: desprezivel. Arredondar para a centena mais
proxima representa erro relativo de no maximo cerca de 2,6 por cento na menor
unidade, de 1.910 ha, e muito menos em todas as outras. Nenhuma conclusao deve
mudar. O ganho e de rigor e de rastreabilidade, nao de resultado.

--------------------------------------------------------------------------------
COMO USAR

1. Abra o ArcGIS Pro 3.1, va em Analysis > Python > Python Window, ou use o
   Notebook do projeto.
2. Ajuste as tres variaveis do bloco CONFIGURACAO abaixo, se necessario.
   Por padrao o script ja aponta para as pastas deste projeto.
3. Cole o script inteiro e execute.
4. Ele grava `Dados/Processados/AOH_area_full_precision.csv` e imprime uma
   comparacao unidade a unidade contra os valores atuais.
5. Depois disso, rode `run_all.R` novamente. O script 01 detecta o arquivo
   sozinho e passa a usar a area de precisao total; se o arquivo nao existir,
   ele segue com os valores atuais e avisa.

O script nao altera nenhum shapefile. Ele trabalha sobre copias em memoria.
--------------------------------------------------------------------------------
"""

import os
import csv
import arcpy

# ==============================================================================
# CONFIGURACAO
# ==============================================================================
PROJECT_ROOT = r"D:\Duda_Nacif_TCC"

# Pasta com os nove shapefiles de unidades amostrais, um por taxon, ja separados
# em partes disjuntas (singlepart). Sao os mesmos que o script 01 do R le para
# extrair os centroides.
SHP_DIR = os.path.join(PROJECT_ROOT, "08_Dados_Especies",
                       "Dados_geo_especies", "Sp_data_singlepart")

OUT_CSV = os.path.join(PROJECT_ROOT, "Dados", "Processados",
                       "AOH_area_full_precision.csv")

# Arquivo atual, usado apenas para a comparacao impressa no final. Se nao
# existir, o script roda do mesmo jeito e pula a comparacao.
CURRENT_CSV = os.path.join(PROJECT_ROOT, "Dados", "Processados",
                           "Data_Raw_FINAL.csv")

# South America Albers Equal Area Conic, o mesmo sistema declarado na Secao 2.4
# do manuscrito. Em uma projecao equivalente a area planar e a medida correta.
ALBERS_WKID = 102033

arcpy.env.overwriteOutput = True

# ==============================================================================
# 1. LOCALIZAR OS SHAPEFILES
# ==============================================================================
if not os.path.isdir(SHP_DIR):
    raise RuntimeError("Pasta de shapefiles nao encontrada: {}".format(SHP_DIR))

shapefiles = sorted(
    os.path.join(SHP_DIR, f) for f in os.listdir(SHP_DIR)
    if f.lower().endswith(".shp")
)
if not shapefiles:
    raise RuntimeError("Nenhum .shp em: {}".format(SHP_DIR))

print("Shapefiles encontrados: {}".format(len(shapefiles)))
for s in shapefiles:
    print("   {}".format(os.path.basename(s)))

albers = arcpy.SpatialReference(ALBERS_WKID)
print("\nSistema de destino: {}".format(albers.name))


# ==============================================================================
# 2. FUNCOES AUXILIARES
# ==============================================================================
def find_field(fields, patterns):
    """Devolve o primeiro campo cujo nome bate com um dos padroes, sem
    diferenciar maiusculas. Devolve None se nenhum bater."""
    lower = {f.name.lower(): f.name for f in fields}
    for p in patterns:
        if p in lower:
            return lower[p]
    for p in patterns:
        for name_l, name in lower.items():
            if p in name_l:
                return name
    return None


def species_from_filename(path):
    """Nome do taxon a partir do nome do arquivo, tirando o sufixo de projecao.
    So e usado quando o shapefile nao tem campo de especie."""
    base = os.path.splitext(os.path.basename(path))[0]
    for suffix in ("_utm", "_albers", "_wgs", "_UTM", "_ALBERS", "_WGS"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


# ==============================================================================
# 3. CALCULAR AS AREAS
# ==============================================================================
rows = []
for shp in shapefiles:
    name = os.path.basename(shp)
    fields = arcpy.ListFields(shp)

    sp_field = find_field(fields, ["sp_id", "species"])
    ua_field = find_field(fields, ["orig_fid", "origfid", "orig_id"])
    if ua_field is None:
        raise RuntimeError(
            "Campo ORIG_FID nao encontrado em {}. Campos disponiveis: {}"
            .format(name, ", ".join(f.name for f in fields))
        )

    # Copia em memoria, reprojetada. O shapefile original nao e tocado.
    tmp = "memory\\aoh_{}".format(
        os.path.splitext(name)[0].replace("-", "_").replace(".", "_")
    )
    src_sr = arcpy.Describe(shp).spatialReference
    if src_sr.factoryCode != ALBERS_WKID:
        arcpy.management.Project(shp, tmp, albers)
    else:
        arcpy.management.CopyFeatures(shp, tmp)

    # Duas medidas. Em uma projecao equivalente as duas praticamente coincidem;
    # a diferenca entre elas e uma verificacao util do sistema de coordenadas.
    arcpy.management.AddField(tmp, "AREA_PLN", "DOUBLE")
    arcpy.management.AddField(tmp, "AREA_GEO", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(
        tmp,
        [["AREA_PLN", "AREA"], ["AREA_GEO", "AREA_GEODESIC"]],
        area_unit="HECTARES",
        coordinate_system=albers,
    )

    read_fields = [ua_field, "AREA_PLN", "AREA_GEO"]
    if sp_field:
        read_fields.insert(0, sp_field)

    n = 0
    with arcpy.da.SearchCursor(tmp, read_fields) as cur:
        for rec in cur:
            if sp_field:
                sp, ua, a_pln, a_geo = rec
            else:
                sp = species_from_filename(shp)
                ua, a_pln, a_geo = rec
            rows.append({
                "SPECIES": str(sp).strip(),
                "UA_ID": int(ua),
                "area_planar_ha": float(a_pln),
                "area_geodesic_ha": float(a_geo),
                "source_shapefile": name,
            })
            n += 1
    arcpy.management.Delete(tmp)
    print("   {}: {} unidade(s)".format(name, n))

print("\nTotal de unidades: {}".format(len(rows)))
if len(rows) != 67:
    print("ATENCAO: o desenho tem 67 unidades e este script leu {}. "
          "Confira se a pasta contem exatamente os nove shapefiles "
          "singlepart usados na analise.".format(len(rows)))


# ==============================================================================
# 4. GRAVAR O CSV COM PRECISAO TOTAL
# ==============================================================================
# repr() do float grava todas as casas significativas. Nao arredonde aqui: o
# objetivo do arquivo e justamente nao perder precisao.
out_dir = os.path.dirname(OUT_CSV)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

rows.sort(key=lambda r: (r["SPECIES"], r["UA_ID"]))
with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["SPECIES", "UA_ID", "area_planar_ha",
                "area_geodesic_ha", "source_shapefile"])
    for r in rows:
        w.writerow([r["SPECIES"], r["UA_ID"],
                    repr(r["area_planar_ha"]), repr(r["area_geodesic_ha"]),
                    r["source_shapefile"]])
print("\nGravado: {}".format(OUT_CSV))


# ==============================================================================
# 5. COMPARACAO COM OS VALORES ATUAIS
# ==============================================================================
if not os.path.isfile(CURRENT_CSV):
    print("\n{} nao encontrado; comparacao pulada.".format(CURRENT_CSV))
else:
    current = {}
    with open(CURRENT_CSV, "r", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        area_col = None
        for cand in ("AOH_unit_area_ha", "AOH_ha", "UA_Area_ha", "EHA_ha"):
            if rd.fieldnames and cand in rd.fieldnames:
                area_col = cand
                break
        if area_col is None:
            print("\nColuna de area nao encontrada em {}. Campos: {}"
                  .format(os.path.basename(CURRENT_CSV), rd.fieldnames))
        else:
            for rec in rd:
                try:
                    key = (str(rec["SPECIES"]).strip(), int(float(rec["UA_ID"])))
                    current[key] = float(rec[area_col])
                except (KeyError, ValueError, TypeError):
                    continue

    if current:
        print("\n{:<16} {:>4} {:>14} {:>16} {:>10}".format(
            "SPECIES", "UA", "atual (ha)", "precisao total", "dif %"))
        print("-" * 66)
        diffs, unmatched = [], 0
        for r in rows:
            key = (r["SPECIES"], r["UA_ID"])
            old = current.get(key)
            if old is None:
                unmatched += 1
                continue
            new = r["area_planar_ha"]
            pct = 100.0 * (new - old) / old if old else float("nan")
            diffs.append(abs(pct))
            print("{:<16} {:>4d} {:>14.1f} {:>16.4f} {:>9.3f}%".format(
                r["SPECIES"][:16], r["UA_ID"], old, new, pct))
        if diffs:
            print("-" * 66)
            print("diferenca absoluta: mediana {:.3f}%  maxima {:.3f}%".format(
                sorted(diffs)[len(diffs) // 2], max(diffs)))
        if unmatched:
            print("{} unidade(s) sem correspondencia por (SPECIES, UA_ID). "
                  "Isso costuma ser diferenca de grafia do taxon; o script 01 "
                  "do R normaliza os nomes ao ler o arquivo.".format(unmatched))

print("\nPronto. Rode run_all.R em seguida.")
