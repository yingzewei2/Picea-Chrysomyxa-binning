# -*- coding: utf-8 -*-
"""
转录组宿主/病原分箱脚本
用于 Picea schrenkiana x Chrysomyxa deformans 双物种混合转录组

输入：欧易交付的 counts_anno.xls（制表符分隔文本，counts + 注释合并表）
输出：宿主 counts 矩阵、病原 counts 矩阵、分箱明细表、待人工判断属名表、统计摘要

判据：
  A  无症状组（H）各样本 counts 全为 0  -> 判为病原
  B  NR 最佳匹配的分类学归属（植物 / 真菌）
  两判据冲突时以注释（判据 B）为准，并在明细表中标记
  无 NR 注释者仅能用判据 A，标记为低置信
"""

import sys

# ---------- 5. Windows 编码适配 ----------
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------- 4. 依赖库检查 ----------
try:
    import pandas as pd
except ModuleNotFoundError:
    print("=" * 60)
    print("缺少依赖库：pandas")
    print("请在命令行运行以下命令安装后重试：")
    print()
    print("    pip install pandas")
    print()
    print("=" * 60)
    try:
        input("按回车键退出...")
    except EOFError:
        pass
    sys.exit(1)

import re
from pathlib import Path
from collections import Counter

# ============================================================
# 分类名录
#   取自本项目 NR 注释中出现的全部属名（前 200 属覆盖 98.4%）
#   未收录的属会写入 待人工判断属名.txt，请人工核对后补入
# ============================================================

PLANT_GENERA = {
    # 裸子植物
    "Picea", "Pinus", "Larix", "Abies", "Taxus", "Ginkgo", "Cycas", "Cunninghamia",
    "Cryptomeria", "Pseudotsuga", "Cephalotaxus", "Juniperus", "Thuja", "Taxodium",
    "Hesperocyparis", "Microbiota", "Gnetum", "Ephedra", "Sequoia", "Metasequoia",
    "Araucaria", "Podocarpus", "Cedrus", "Tsuga", "Sciadopitys", "Welwitschia",
    # 蕨类、苔藓、藻类
    "Adiantum", "Ceratopteris", "Selaginella", "Marchantia", "Sphagnum", "Ceratodon",
    "Physcomitrium", "Physcomitrella", "Azolla", "Salvinia", "Anthoceros", "Chara",
    "Chlamydomonas", "Klebsormidium", "Volvox", "Ostreococcus", "Micromonas",
    # 基部被子植物
    "Amborella", "Nymphaea", "Cinnamomum", "Aristolochia", "Tetracentron", "Kingdonia",
    "Nelumbo", "Macleaya", "Papaver", "Coptis", "Thalictrum", "Aquilegia", "Telopea",
    "Macadamia", "Trochodendron", "Buxus", "Liriodendron", "Magnolia", "Piper",
    # 真双子叶、单子叶
    "Vitis", "Gossypium", "Prunus", "Quercus", "Populus", "Salix", "Nyssa", "Actinidia",
    "Trifolium", "Citrus", "Lupinus", "Solanum", "Camellia", "Carya", "Syzygium",
    "Rhododendron", "Tanacetum", "Shorea", "Juglans", "Durio", "Theobroma", "Nicotiana",
    "Vigna", "Brassica", "Castanea", "Corchorus", "Cichorium", "Rhamnella", "Punica",
    "Bauhinia", "Senna", "Vaccinium", "Hevea", "Parasponia", "Morus", "Olea", "Datura",
    "Corymbia", "Pyrus", "Sesamum", "Cucurbita", "Jatropha", "Pistacia", "Arachis",
    "Cicer", "Beta", "Arabidopsis", "Malus", "Cajanus", "Ziziphus", "Cucumis", "Glycine",
    "Paulownia", "Smallanthus", "Trema", "Carpinus", "Coffea", "Melastoma", "Mangifera",
    "Phtheirospermum", "Catharanthus", "Ipomoea", "Ricinus", "Handroanthus", "Mikania",
    "Carica", "Cannabis", "Impatiens", "Striga", "Morella", "Perilla", "Chenopodium",
    "Fragaria", "Artemisia", "Prosopis", "Abrus", "Acer", "Xanthoceras", "Pisum",
    "Medicago", "Lotus", "Phaseolus", "Helianthus", "Lactuca", "Daucus", "Capsicum",
    "Eucalyptus", "Betula", "Fagus", "Ulmus", "Fraxinus", "Rosa", "Rubus", "Pelargonium",
    "Elaeis", "Phoenix", "Musa", "Dioscorea", "Ananas", "Hordeum", "Cocos", "Colocasia",
    "Dendrobium", "Zingiber", "Ensete", "Asparagus", "Triticum", "Zea", "Oryza",
    "Panicum", "Phalaenopsis", "Apostasia", "Vanilla", "Digitaria", "Spirodela",
    "Setaria", "Sorghum", "Brachypodium", "Aegilops", "Lolium", "Amaranthus",
    # 补充（首轮运行后由 待人工判断属名.txt 核对补入）
    "Momordica", "Cuscuta", "Buddleja", "Mucuna", "Hibiscus", "Herrania", "Spatholobus",
    "Cupressus", "Callitropsis", "Lithospermum", "Camelina", "Arabis", "Dorcoceras",
    "Tripterygium", "Salvia", "Eleusine", "Cynara", "Erigeron", "Psidium", "Erythranthe",
    "Manihot", "Miscanthus", "Calocedrus", "Spinacia", "Carex", "Cephalotus", "Rhodamnia",
    "Genlisea", "Capsella", "Zizania", "Eragrostis", "Tarenaya", "Amentotaxus", "Arctium",
    "Benincasa", "Thlaspi", "Taiwania", "Heliosperma", "Microthlaspi", "Nothotsuga",
    "Raphanus", "Vernicia", "Zostera", "Conocephalum", "Sesbania", "Aethionema",
    # 第二轮补充（由 待人工判断属名.txt 核对补入）
    "Consolida", "Diospyros", "Nageia", "Agave", "Funaria", "Lagenaria", "Sorbus",
    "Oenothera", "Pohlia", "Paeonia", "Chimonanthus", "Hippophae", "Narcissus",
    "Takakia", "Chamaecyparis", "Itatiella", "Ruscus", "Aquilaria", "Akebia", "Annona",
    "Chrysanthemum", "Canella", "Panax", "Chelidonium", "Hosta", "Petunia", "Apium",
    "Dichanthelium", "Grimmia", "Apopellia", "Auxenochlorella", "Ulva", "Tetrabaena",
}

FUNGAL_GENERA = {
    # 锈菌（柄锈菌目）—— 与本研究病原同目，是病原侧的主要比对目标
    "Cronartium", "Puccinia", "Austropuccinia", "Chrysomyxa", "Melampsora", "Uromyces",
    "Phakopsora", "Gymnosporangium", "Hemileia", "Coleosporium", "Pucciniastrum",
    "Milesina", "Thekopsora", "Uredinopsis", "Endocronartium", "Peridermium",
    "Sphaerophragmium", "Tranzschelia", "Physopella", "Rossmanomyces",
    # 座囊菌纲等其他子囊菌
    "Pyrenochaeta", "Alternaria", "Stagonospora", "Plenodomus", "Sclerotinia",
    "Phaeosphaeria", "Pyrenophora", "Botrytis", "Corynespora", "Bipolaris",
    "Lentithecium", "Melanomma", "Bimuria", "Rutstroemia", "Karstenula", "Zopfia",
    "Massarina", "Parastagonospora", "Lophiostoma", "Rachicladosporium", "Pleomassaria",
    "Polychaeton", "Dothistroma", "Fusarium", "Aureobasidium", "Amniculicola",
    "Sporormia", "Cercospora", "Cadophora", "Delitschia", "Colletotrichum",
    "Hyaloscypha", "Pseudogymnoascus", "Saccharata", "Oidiodendron", "Hortaea",
    "Penicillium", "Aspergillus", "Lepidopterella", "Neofusicoccum", "Bisporella",
    "Cladosporium", "Verticillium", "Magnaporthe", "Trichoderma", "Neurospora",
    "Septoria", "Zymoseptoria", "Venturia", "Diaporthe", "Phomopsis", "Epicoccum",
    "Didymella", "Ascochyta", "Leptosphaeria", "Ophiostoma", "Grosmannia",
    "Heterobasidion", "Armillaria", "Rhizoctonia", "Sclerotium",
    # 担子菌（非锈菌）
    "Cryptococcus", "Saitozyma", "Kwoniella", "Rhodotorula", "Ceratobasidium",
    "Tulasnella", "Serendipita", "Termitomyces", "Peniophora", "Mycena", "Amanita",
    "Calocera", "Plicaturopsis", "Microbotryum", "Malassezia", "Ustilago", "Tilletia",
    "Sporobolomyces", "Trametes", "Ganoderma", "Schizophyllum", "Postia", "Laccaria",
    # 科/目级名称（NR 中偶见）
    "Phaeosphaeriaceae", "Tremellales", "Helotiales", "Xylariaceae", "Pleosporaceae",
    "Dothideomycetes", "Sordariomycetes", "Leotiomycetes", "Agaricomycetes",
    "Pucciniales", "Pucciniomycetes", "Ascomycota", "Basidiomycota", "Fungi",
    # 补充（首轮运行后由 待人工判断属名.txt 核对补入）
    "Leptodontidium", "Auricularia", "Dendrothele", "Myriangium", "Macrophomina",
    "Galerina", "Cylindrobasidium", "Botryobasidium", "Jaapia", "Moniliophthora",
    "Fibularhizoctonia", "Coprinopsis", "Xylariales", "Patellaria", "Rhizophagus",
    "Mixia", "Podila", "Aulographum", "Acidomyces", "Neopestalotiopsis", "Coniochaeta",
    "Diaporthaceae", "Rhexocercosporidium", "Hypoxylon", "Dentipellis", "Xanthoria",
    "Flagelloscypha", "Xylogone", "Nemania", "Xylaria", "Neolentinus", "Sphaerobolus",
    "Agaricus", "Athelia", "Piedraia", "Clavulina", "Gymnopus", "Monosporascus",
    "Eutypa", "Chalara", "Cenococcum", "Stachybotrys", "Paxillus", "Hydnum", "Psilocybe",
    "Macrolepiota", "Tricholoma", "Glarea", "Volvariella", "Hydnomerulius", "Taphrina",
    "Pleurotus", "Gelatoporia", "Coccidioides", "Hebeloma", "Kalmanozyma", "Emergomyces",
    "Lipomyces", "Claviceps", "Suillus", "Auriculariales", "Hypholoma", "Melanopsichium",
    "Vararia", "Talaromyces", "Fistulina", "Pyricularia", "Panus", "Cordyceps",
    "Mariannaea", "Trichosporon", "Piloderma", "Chaetothyriales", "Ilyonectria",
    "Rhizopus", "Mucor", "Mortierella", "Umbelopsis", "Saccharomyces", "Candida",
    # 第二轮补充（由 待人工判断属名.txt 核对补入）
    "Marasmius", "Blumeria", "Wolfiporia", "Morchella", "Spizellomyces", "Epichloe",
    "Wilcoxina", "Apophysomyces", "Paecilomyces", "Neolecta", "Akanthomyces", "Rozella",
    "Chaetosphaeronema", "Claroideoglomus", "Clonostachys", "Russula", "Histoplasma",
    "Basidiobolus", "Nothophoma", "Xanthomendoza", "Pestalotiopsis", "Phlebopus",
    "Paracoccidioides", "Exophiala", "Pyronema", "Rhizopogon", "Tephrocybe", "Chaetomium",
    "Annulohypoxylon", "Yarrowia", "Scleroderma", "Serpula", "Endogone", "Beauveria",
    "Daedalea", "Ophiocordyceps", "Daldinia", "Dactylellina", "Magnaporthiopsis",
    "Pisolithus", "Cubamyces", "Thozetella", "Helicocarpus", "Acarospora", "Trichophyton",
    "Phlebiopsis", "Terfezia", "Saitoella", "Lichtheimia", "fungal",
}

# 既非植物也非真菌（卵菌、细菌、动物等），单独归类，不计入宿主
OTHER_GENERA = {
    "Phytophthora", "Pythium", "Plasmopara", "Albugo", "Peronospora", "Saprolegnia",
    "Escherichia", "Pseudomonas", "Bacillus", "Agrobacterium", "Xanthomonas",
    "Streptomyces", "Rhizobium", "Homo", "Mus", "Drosophila", "Caenorhabditis",
    "Aphanomyces", "Bremia", "Hyaloperonospora", "Thalassiosira", "Phaeodactylum",
    # 第二轮补充：细菌、藻类、原生生物、微孢子虫，以及 NR 中的非物种写法
    "Clostridium", "Symbiodinium", "Porphyridium", "Pelagomonas", "Fistulifera",
    "Lingulodinium", "Tribonema", "Amoeboaphelidium", "Nematocida", "synthetic",
}

# 本研究关注的主线基因关键词，用于分箱后自查
FOCUS_KEYWORDS = [
    # 赤霉素
    "gibberellin", "GA20ox", "GA2ox", "GA3ox", "ent-kaurene", "GID1",
    # 二萜树脂酸
    "diterpene synthase", "levopimaradiene", "abietadiene", "isopimaradiene",
    "CYP720B", "abietadienol",
    # 前体通路
    "1-deoxy-D-xylulose", "DXS", "DXR",
    "hydroxymethylglutaryl", "HMG-CoA", "mevalonate",
    # 类黄酮上游
    "chalcone synthase", "chalcone isomerase", "flavanone 3", "flavonol synthase",
    "phenylalanine ammonia", "cinnamate 4-hydroxylase", "4-coumarate",
    # 缩合单宁 / 原花青素分支（新增）
    "dihydroflavonol 4-reductase", "leucoanthocyanidin reductase",
    "anthocyanidin reductase", "anthocyanidin synthase",
    "leucoanthocyanidin dioxygenase", "flavonoid 3'-hydroxylase",
    # 氧化偶联酶（新增，对应脱氢双儿茶素 A）
    "peroxidase", "polyphenol oxidase", "laccase", "catechol oxidase",
    # 芪类（新增）
    "stilbene synthase", "resveratrol synthase", "pinosylvin",
    # 光合与叶绿素
    "chlorophyll a-b binding", "stay-green", "pheophytinase", "chlorophyllase",
    "ribulose bisphosphate",
    # 表皮屏障（新增）
    "CYP86", "CYP703", "fatty acyl-CoA reductase", "wax synthase",
    "3-ketoacyl-CoA synthase", "eceriferum",
]

# ---------- 判据 C：GC 含量（可选，需 Unigene.fa） ----------
# 本项目实测：植物 Unigene 中位 GC 42.2%，真菌 49.4%，重叠较多。
# 单一阈值的总体准确率仅约 83%，故不作为独立判据，
# 只在两端的高置信区间内使用（实测该区间内准确率 >=92%）。
GC_PLANT_MAX = 42.0      # GC 低于此值：实测仅 8% 为真菌 -> 支持判为宿主
GC_FUNGUS_MIN = 52.0     # GC 高于此值：实测 93% 为真菌 -> 支持判为病原
GC_MIN_LENGTH = 300      # 序列过短时 GC 不稳定，低于此长度不使用判据 C

H_PREFIX = "H-"          # 无症状组样本列的前缀
SAMPLE_PREFIXES = ("H-", "S1-", "S2-")


def find_input_file(script_dir: Path) -> Path:
    """查找 counts_anno.xls。依次尝试：脚本目录 -> 当前目录 -> 交互输入。
    交互输入既接受文件路径，也接受文件夹路径（会在其中及子目录中查找）。"""
    names = ["counts_anno.xls", "counts_anno.txt", "counts_anno.tsv", "counts_anno.csv"]

    def look_in(folder: Path):
        """在文件夹本身及其子目录中查找目标文件。"""
        for n in names:
            p = folder / n
            if p.is_file():
                return p
        for n in names:
            try:
                hits = sorted(folder.rglob(n))
            except (PermissionError, OSError):
                hits = []
            if hits:
                return hits[0]
        return None

    for d in (script_dir, Path.cwd()):
        hit = look_in(d)
        if hit:
            return hit

    print("未在脚本所在目录找到 counts_anno.xls。")
    print()
    print("请粘贴 counts_anno.xls 的路径，或它所在文件夹的路径，然后按回车。")
    print("（可以直接把文件或文件夹拖到本窗口里，路径会自动填入）")
    print("直接按回车则退出。")
    print()

    for attempt in range(5):
        raw = input("路径> ").strip().strip('"').strip("'").rstrip("\\")
        if not raw:
            raise SystemExit("已取消。")
        p = Path(raw)
        if p.is_file():
            if p.name.lower() not in [n.lower() for n in names]:
                print(f"  注意：文件名是 {p.name}，不是 counts_anno.xls，仍按该文件读取。")
            return p
        if p.is_dir():
            hit = look_in(p)
            if hit:
                print(f"  在该文件夹中找到：{hit}")
                return hit
            print(f"  [未找到] 文件夹存在，但其中没有 counts_anno.xls：{p}")
            print("  该文件夹下的 .xls / .txt 文件有：")
            found = False
            for f in sorted(p.glob("*.xls")) + sorted(p.glob("*.txt")) + sorted(p.glob("*.tsv")):
                print(f"    {f.name}")
                found = True
            if not found:
                print("    （没有）")
            print("  请改为直接粘贴目标文件的完整路径。")
        else:
            print(f"  [无效] 该路径不存在：{p}")
            print("  提示：路径中若含空格或中文，直接把文件拖进窗口最稳妥。")
        print()

    raise FileNotFoundError("多次尝试均未找到有效文件，已退出。")


def read_fasta_gc(path: Path):
    """读取 FASTA，返回 {序列ID: (GC%, 长度)}。逐行处理，不占大内存。"""
    gc, ident, n_gc, n_len = {}, None, 0, 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if ident is not None:
                    gc[ident] = (100.0 * n_gc / n_len if n_len else None, n_len)
                ident = line[1:].split()[0]
                n_gc = n_len = 0
            else:
                seq = line.strip().upper()
                n_len += len(seq)
                n_gc += seq.count("G") + seq.count("C")
        if ident is not None:
            gc[ident] = (100.0 * n_gc / n_len if n_len else None, n_len)
    return gc


def find_fasta(script_dir: Path, near: Path):
    """在脚本目录、输入文件同级及其上级目录中查找 Unigene.fa（可选）。"""
    names = ["Unigene.fa", "Unigene.fasta", "unigene.fa"]
    roots = [script_dir, near.parent, near.parent.parent,
             near.parent.parent / "1.2.De_novo", Path.cwd()]
    for d in roots:
        try:
            for n in names:
                p = d / n
                if p.is_file():
                    return p
        except OSError:
            continue
    return None


def parse_species(nr_annotation: str) -> str:
    """从 NR 注释中取方括号内的物种名，取最后一个方括号。"""
    if not isinstance(nr_annotation, str) or nr_annotation.strip() in ("", "--", "-"):
        return ""
    hits = re.findall(r"\[([^\[\]]+)\]", nr_annotation)
    return hits[-1].strip() if hits else ""


def classify_genus(species: str):
    """返回 (属名, 分类)。分类为 plant / fungus / other / unknown / none。"""
    if not species:
        return "", "none"
    genus = species.split()[0]
    if genus in PLANT_GENERA:
        return genus, "plant"
    if genus in FUNGAL_GENERA:
        return genus, "fungus"
    if genus in OTHER_GENERA:
        return genus, "other"
    return genus, "unknown"


def main():
    script_dir = Path(__file__).resolve().parent
    infile = find_input_file(script_dir)
    outdir = script_dir / "分箱结果"
    outdir.mkdir(exist_ok=True)

    print(f"读取：{infile}")
    print("文件较大（约 35 MB），请稍候…")
    df = pd.read_csv(infile, sep="\t", dtype={"id": str}, low_memory=False)

    # 去掉尾部空列（原表行末有多余制表符）
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed:")]]

    if "id" not in df.columns:
        raise ValueError("表中未找到 id 列，请确认输入文件是否为 counts_anno.xls")
    if "NR annotation" not in df.columns:
        raise ValueError("表中未找到 'NR annotation' 列，请确认输入文件")

    sample_cols = [c for c in df.columns if str(c).startswith(SAMPLE_PREFIXES)]
    h_cols = [c for c in sample_cols if str(c).startswith(H_PREFIX)]
    if not h_cols:
        raise ValueError("未找到无症状组（H-）样本列，判据 A 无法执行")

    print(f"Unigene 数：{len(df)}")
    print(f"样本列（{len(sample_cols)} 个）：{', '.join(sample_cols)}")
    print(f"无症状组（判据 A 使用）：{', '.join(h_cols)}")
    print()

    # ---------- 判据 A ----------
    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df[sample_cols] = counts
    h_all_zero = (counts[h_cols] == 0).all(axis=1)
    df["判据A_H组全零"] = h_all_zero
    df["判据A_结论"] = h_all_zero.map({True: "pathogen", False: "host"})

    # ---------- 判据 C：GC 含量（可选） ----------
    fa = find_fasta(script_dir, infile)
    use_gc = False
    if fa is not None:
        print(f"找到组装序列：{fa}")
        print("正在计算 GC 含量…")
        gcmap = read_fasta_gc(fa)
        df["GC含量"] = df["id"].map(lambda x: (gcmap.get(x) or (None, 0))[0])
        df["序列长度"] = df["id"].map(lambda x: (gcmap.get(x) or (None, 0))[1])
        matched = int(df["GC含量"].notna().sum())
        print(f"  匹配到 GC 的 Unigene：{matched} / {len(df)}")
        if matched > 0.5 * len(df):
            use_gc = True
        else:
            print("  [警告] 匹配率过低，判据 C 不启用（请确认 fa 与 counts 表来自同一次组装）")
    else:
        print("未找到 Unigene.fa，跳过判据 C（GC 含量）。")
        print("  如需启用，把 Unigene.fa 放到本脚本或 counts_anno.xls 所在目录。")
        df["GC含量"] = None
        df["序列长度"] = 0
    print()

    if use_gc:
        def gc_call(row):
            g, L = row["GC含量"], row["序列长度"]
            if g is None or pd.isna(g) or L < GC_MIN_LENGTH:
                return ""
            if g < GC_PLANT_MAX:
                return "host"
            if g > GC_FUNGUS_MIN:
                return "pathogen"
            return ""
        df["判据C_结论"] = df.apply(gc_call, axis=1)
    else:
        df["判据C_结论"] = ""

    # ---------- 判据 B ----------
    species = df["NR annotation"].map(parse_species)
    parsed = species.map(classify_genus)
    df["NR物种"] = species
    df["NR属名"] = parsed.map(lambda x: x[0])
    df["分类归属"] = parsed.map(lambda x: x[1])
    df["判据B_结论"] = df["分类归属"].map(
        {"plant": "host", "fungus": "pathogen", "other": "other",
         "unknown": "", "none": ""}
    )

    # ---------- 合并判据 ----------
    final, conf, conflict = [], [], []
    for b, a, c, cls in zip(df["判据B_结论"], df["判据A_结论"],
                            df["判据C_结论"], df["分类归属"]):
        if b in ("host", "pathogen"):
            # 有 NR 注释：以注释为准（判据 B 优先级最高）
            final.append(b)
            conf.append("high")
            conflict.append(b != a)
        elif b == "other":
            final.append("other")
            conf.append("high")
            conflict.append(False)
        elif c in ("host", "pathogen"):
            # 无注释，但 GC 落在高置信区间：A 与 C 一致则中置信，冲突则以 A 为准
            if c == a:
                final.append(a)
                conf.append("medium")
                conflict.append(False)
            else:
                final.append(a)
                conf.append("low")
                conflict.append(True)
        else:
            final.append(a)
            conf.append("low")          # 仅判据 A，无任何佐证
            conflict.append(False)
    df["最终归属"] = final
    df["置信度"] = conf
    df["判据冲突"] = conflict

    # ---------- 输出 ----------
    host = df[df["最终归属"] == "host"]
    patho = df[df["最终归属"] == "pathogen"]
    other = df[df["最终归属"] == "other"]

    host_mat = host[["id"] + sample_cols].copy()
    patho_mat = patho[["id"] + sample_cols].copy()

    host_mat.to_csv(outdir / "host_counts.tsv", sep="\t", index=False)
    patho_mat.to_csv(outdir / "pathogen_counts.tsv", sep="\t", index=False)

    # DESeq2 需要整数，另出取整版
    host_int = host_mat.copy()
    host_int[sample_cols] = host_int[sample_cols].round().astype("int64")
    host_int.to_csv(outdir / "host_counts_integer_for_DESeq2.tsv", sep="\t", index=False)

    # 宿主子集内重新归一化的 CPM（供快速查看；正式差异分析仍交 DESeq2）
    cpm = host_mat.copy()
    lib = cpm[sample_cols].sum(axis=0)
    cpm[sample_cols] = cpm[sample_cols].div(lib, axis=1) * 1e6
    cpm.to_csv(outdir / "host_CPM_renormalised.tsv", sep="\t", index=False)

    detail_cols = ["id", "NR物种", "NR属名", "分类归属", "GC含量", "序列长度",
                   "判据A_H组全零", "判据A_结论", "判据B_结论", "判据C_结论",
                   "最终归属", "置信度", "判据冲突"]
    detail_cols = [c for c in detail_cols if c in df.columns]
    df[detail_cols].to_csv(outdir / "分箱明细表.tsv", sep="\t", index=False)

    # 待人工判断的属名
    unk = Counter(df.loc[df["分类归属"] == "unknown", "NR属名"])
    with open(outdir / "待人工判断属名.txt", "w", encoding="utf-8") as f:
        f.write("以下属名未收录在脚本的分类名录中，已按判据 A 处理并标记为低置信。\n")
        f.write("请人工核对后，把它们补进脚本顶部的 PLANT_GENERA / FUNGAL_GENERA / OTHER_GENERA，重跑一次。\n")
        f.write("（按 Unigene 数从多到少排列）\n\n")
        for g, n in unk.most_common():
            f.write(f"{g}\t{n}\n")

    # ---------- 主线基因自查 ----------
    ann = (df["NR annotation"].fillna("") + " | " +
           df.get("Swissprot annotation", pd.Series([""] * len(df))).fillna("") + " | " +
           df.get("KEGG description", pd.Series([""] * len(df))).fillna(""))
    focus_rows = []
    for kw in FOCUS_KEYWORDS:
        hit = df[ann.str.contains(re.escape(kw), case=False, na=False)]
        if len(hit):
            vc = hit["最终归属"].value_counts().to_dict()
            focus_rows.append({
                "关键词": kw, "命中数": len(hit),
                "宿主": vc.get("host", 0), "病原": vc.get("pathogen", 0),
                "其他": vc.get("other", 0),
            })
    focus = pd.DataFrame(focus_rows)
    if len(focus):
        focus.to_csv(outdir / "主线基因归属自查.tsv", sep="\t", index=False)

    # 逐条明细：每个命中的 Unigene 属于宿主还是病原，便于人工核对
    mask = pd.Series(False, index=df.index)
    for kw in FOCUS_KEYWORDS:
        mask |= ann.str.contains(re.escape(kw), case=False, na=False)
    hit_cols = ["id", "NR物种", "最终归属", "置信度",
                "NR annotation", "Swissprot annotation"]
    hit_cols = [c for c in hit_cols if c in df.columns]
    df.loc[mask, hit_cols + sample_cols].to_csv(
        outdir / "主线基因逐条明细.tsv", sep="\t", index=False)

    # ---------- 摘要 ----------
    rust = {"Cronartium", "Puccinia", "Austropuccinia", "Chrysomyxa", "Melampsora",
            "Uromyces", "Phakopsora", "Gymnosporangium", "Hemileia", "Coleosporium",
            "Pucciniastrum", "Milesina", "Thekopsora", "Uredinopsis", "Rossmanomyces"}
    n_rust = int(df["NR属名"].isin(rust).sum())
    n_fungus = int((df["分类归属"] == "fungus").sum())

    lines = []
    A = lines.append
    A("=" * 62)
    A("宿主 / 病原 分箱结果摘要")
    A("=" * 62)
    A(f"输入文件              : {infile.name}")
    A(f"Unigene 总数          : {len(df)}")
    A(f"样本数                : {len(sample_cols)}（无症状组 {len(h_cols)} 个）")
    A("")
    A("— NR 注释情况 —")
    A(f"  有 NR 物种信息      : {int((species != '').sum())} ({100*(species!='').mean():.1f}%)")
    A(f"  无 NR 注释          : {int((species == '').sum())} ({100*(species=='').mean():.1f}%)")
    A(f"  注释为植物          : {int((df['分类归属']=='plant').sum())}")
    A(f"  注释为真菌          : {n_fungus}（其中锈菌类 {n_rust}，非锈菌真菌 {n_fungus-n_rust}）")
    A(f"  注释为其他生物      : {int((df['分类归属']=='other').sum())}")
    A(f"  属名未收录          : {int((df['分类归属']=='unknown').sum())}  -> 见 待人工判断属名.txt")
    A("")
    A("— 分箱结果 —")
    A(f"  宿主 Unigene        : {len(host)} ({100*len(host)/len(df):.1f}%)")
    A(f"  病原 Unigene        : {len(patho)} ({100*len(patho)/len(df):.1f}%)")
    A(f"  其他生物            : {len(other)}")
    A(f"  判据冲突（以注释为准）: {int(df['判据冲突'].sum())} ({100*df['判据冲突'].mean():.1f}%)")
    A(f"  中置信（判据 A+C 一致）: {int((df['置信度']=='medium').sum())} ({100*(df['置信度']=='medium').mean():.1f}%)")
    A(f"  低置信（仅判据 A）  : {int((df['置信度']=='low').sum())} ({100*(df['置信度']=='low').mean():.1f}%)")
    if use_gc:
        A("")
        A("— 判据 C（GC 含量）—")
        A(f"  已启用。阈值：GC < {GC_PLANT_MAX}% 支持宿主，GC > {GC_FUNGUS_MIN}% 支持病原")
        A(f"  仅用于无 NR 注释者；本项目实测该两端区间内准确率 >= 92%")
        A(f"  由判据 C 提升为中置信 : {int((df['置信度']=='medium').sum())}")
    A("")
    A("— 各样本宿主文库深度（分箱后 counts 之和）—")
    for c in sample_cols:
        A(f"  {c:<8} {host_mat[c].sum():>16,.0f}")
    A("")
    A("— 各样本病原占比（counts 计）—")
    for c in sample_cols:
        tot = counts[c].sum()
        p = patho_mat[c].sum() / tot * 100 if tot else 0
        A(f"  {c:<8} {p:>6.2f}%")
    A("")
    A("输出目录：" + str(outdir))
    A("  host_counts.tsv                     宿主 counts 矩阵（原始小数）")
    A("  host_counts_integer_for_DESeq2.tsv  宿主 counts 取整版，供 DESeq2")
    A("  host_CPM_renormalised.tsv           宿主子集内重新归一化的 CPM")
    A("  pathogen_counts.tsv                 病原 counts 矩阵")
    A("  分箱明细表.tsv                      每条 Unigene 的判据与结论")
    A("  待人工判断属名.txt                  未收录属名，需人工核对后补入脚本")
    A("  主线基因归属自查.tsv                关注基因落在宿主侧还是病原侧")
    A("=" * 62)

    text = "\n".join(lines)
    print(text)
    (outdir / "分箱统计摘要.txt").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
        print()
        print("完成。")
    except Exception as e:
        print()
        print("=" * 60)
        print("运行出错：")
        print(f"  {type(e).__name__}: {e}")
        print()
        import traceback
        traceback.print_exc()
        print("=" * 60)
        print("把上面的报错内容整段复制给我，我来定位问题。")
    try:
        input("按回车键退出...")
    except EOFError:
        pass
