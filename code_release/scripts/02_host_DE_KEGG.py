# -*- coding: utf-8 -*-
"""
宿主子集差异表达分析（配对模型） + 宿主背景 KEGG 富集
输出为欧易云联合分析所需的目录结构与文件格式

解决三件事：
  1. 差异分析在宿主分箱子集内重新归一化，不受病原转录本影响
  2. 采用配对模型 design = ~ tree + stage，用上同株配对设计
  3. KEGG 富集以宿主 Unigene 为背景，而非全组装

输入（放在同一文件夹）：
  host_counts_integer_for_DESeq2.tsv   分箱脚本的输出
  counts_anno.xls                      注释来源（NR、GO、Pfam 等列）
  Unigene.KEGG.gene.anno.xls           KEGG 完整注释（一基因对多通路）
  Unigene.KEGG.classification.xls      KEGG 一级/二级分类
  Unigene.fa                           可选，用于计算宿主内重归一化 FPKM

说明：counts_anno.xls 中每个基因仅保留一条通路，且把 "koA|koB" 当作单一条目，
      会导致富集统计失真。故优先读取 Unigene.KEGG.gene.anno.xls；
      若缺失则回退到 counts_anno.xls 并给出警告。

输出：2.gene/ 目录，结构与欧易报告一致，可直接打包上传
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

_MISSING = []
try:
    import pandas as pd
    import numpy as np
except ModuleNotFoundError as e:
    _MISSING.append(e.name)
try:
    from scipy.stats import hypergeom
except ModuleNotFoundError:
    _MISSING.append("scipy")
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ModuleNotFoundError:
    _MISSING.append("pydeseq2")

if _MISSING:
    print("=" * 62)
    print("缺少依赖库：" + "、".join(sorted(set(_MISSING))))
    print("请在命令行运行以下命令安装后重试：")
    print()
    print("    pip install pandas numpy scipy pydeseq2")
    print()
    print("pydeseq2 会连带安装 anndata、statsmodels、formulaic，可能需要一两分钟。")
    print("=" * 62)
    try:
        input("按回车键退出...")
    except EOFError:
        pass
    sys.exit(1)

import re
import warnings
import logging
from pathlib import Path

warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)

# ============================================================
# 参数
# ============================================================

CONTRASTS = [("S1", "H"), ("S2", "H"), ("S2", "S1")]   # (case, control)
MIN_MEAN_COUNT = 2        # 与欧易一致：counts 均值 <= 2 的 Unigene 剔除
Q_CUT = 0.05
FC_CUT = 2.0              # |FoldChange| >= 2，即 |log2FC| >= 1

ANNO_COLS = [
    "NR id", "NR annotation", "Swissprot id", "Swissprot annotation",
    "KEGG id", "KO", "KEGG gene name", "KEGG description",
    "Pathway", "Pathway definition",
    "KOG id", "KOG defination", "KOG code", "KOG categories", "KOG class",
    "eggNOG abbreviation", "eggNOG id", "eggNOG description",
    "GO annotation", "Pfam IDs", "Pfam Description", "iTAK_type", "Family",
]

def load_kegg(script_dir, counts_anno_path, host_ids):
    """读取 KEGG 注释。优先用完整版（一基因对多通路），否则回退到 counts_anno。
    返回 (pairs, pw_name, pw_cls, source)
      pairs  : [(Unigene, ko), ...]
      pw_name: {ko: 通路名}
      pw_cls : {ko: (一级分类, 二级分类)}
    """
    pairs, pw_name, pw_cls = [], {}, {}

    f = find_file(script_dir, ["Unigene.KEGG.gene.anno.xls",
                               "Unigene_KEGG_gene_anno.xls"], required=False)
    if f is not None:
        g = pd.read_csv(f, sep="\t", low_memory=False)
        g.columns = [str(c).lstrip("#") for c in g.columns]
        for gid, pw, pd_ in zip(g["Gene_id"], g["Pathway"].astype(str),
                                g["Pathway_definition"].astype(str)):
            kos = [x.strip() for x in pw.split("|")]
            nms = [x.strip() for x in pd_.split("|")]
            for i, ko in enumerate(kos):
                if ko.startswith("ko"):
                    pairs.append((gid, ko))
                    if i < len(nms):
                        pw_name[ko] = nms[i]
        source = f"完整 KEGG 注释（{f.name}）"
    else:
        a = pd.read_csv(counts_anno_path, sep="\t", low_memory=False,
                        usecols=["id", "Pathway", "Pathway definition"])
        for gid, pw, pd_ in zip(a["id"], a["Pathway"].astype(str),
                                a["Pathway definition"].astype(str)):
            if pw.startswith("ko") and "|" not in pw:
                pairs.append((gid, pw))
                pw_name[pw] = pd_
        source = ("counts_anno.xls 的简化注释 —— 每基因仅一条通路，"
                  "建议改用 Unigene.KEGG.gene.anno.xls")

    fc = find_file(script_dir, ["Unigene.KEGG.classification.xls",
                                "Unigene_KEGG_classification.xls"], required=False)
    ko_of_gene = {}
    for gid, ko in pairs:
        ko_of_gene.setdefault(gid, set()).add(ko)
    if fc is not None:
        c = pd.read_csv(fc, sep="\t", low_memory=False)
        cnt = {}
        for _, row in c.iterrows():
            l1 = str(row["Classification_level1"]).strip()
            l2 = str(row["Classification_level2"]).strip()
            for gid in str(row["Genes"]).split(";"):
                gid = gid.strip()
                for ko in ko_of_gene.get(gid, ()):
                    cnt.setdefault(ko, {}).setdefault((l1, l2), 0)
                    cnt[ko][(l1, l2)] += 1
        for ko, d in cnt.items():
            pw_cls[ko] = max(d, key=d.get)
    return pairs, pw_name, pw_cls, source


def find_file(script_dir: Path, names, required=True, desc=""):
    for d in (script_dir, Path.cwd()):
        for n in names:
            p = d / n
            if p.is_file():
                return p
        for n in names:
            try:
                hits = sorted(d.rglob(n))
            except OSError:
                hits = []
            if hits:
                return hits[0]
    if not required:
        return None
    print(f"未找到{desc}（{ '、'.join(names) }）。")
    print("请粘贴该文件的完整路径，或它所在文件夹的路径，然后按回车。直接回车退出。")
    for _ in range(5):
        raw = input("路径> ").strip().strip('"').strip("'").rstrip("\\")
        if not raw:
            raise SystemExit("已取消。")
        p = Path(raw)
        if p.is_file():
            return p
        if p.is_dir():
            for n in names:
                q = p / n
                if q.is_file():
                    return q
                hits = sorted(p.rglob(n))
                if hits:
                    return hits[0]
            print(f"  [未找到] 该文件夹中没有 { '、'.join(names) }")
        else:
            print(f"  [无效] 路径不存在：{p}")
    raise FileNotFoundError(f"未找到{desc}")


def read_fasta_lengths(path: Path):
    lens, ident, n = {}, None, 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if ident is not None:
                    lens[ident] = n
                ident = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        if ident is not None:
            lens[ident] = n
    return lens


def enrich(deg_ids, pairs, all_ids, pw_name, pw_cls):
    """超几何富集。pairs: [(Unigene, ko)]，允许一基因对多通路。"""
    bg_pairs = [(g, k) for g, k in pairs if g in all_ids]
    bg_genes = {g for g, _ in bg_pairs}
    fg_genes = bg_genes & set(deg_ids)
    pop_total, list_total = len(bg_genes), len(fg_genes)
    if list_total == 0 or pop_total == 0:
        return pd.DataFrame()
    pop_hits, list_hits = {}, {}
    for g, k in bg_pairs:
        pop_hits.setdefault(k, set()).add(g)
        if g in fg_genes:
            list_hits.setdefault(k, set()).add(g)
    rows = []
    for p, genes in list_hits.items():
        k, n = len(genes), len(pop_hits[p])
        pv = hypergeom.sf(k - 1, pop_total, n, list_total)
        cls = pw_cls.get(p, ("--", "--"))
        rows.append({
            "id": p,
            "Classification_level1": cls[0],
            "Classification_level2": cls[1],
            "Term": pw_name.get(p, "--"),
            "ListHits": k, "ListTotal": list_total,
            "PopHits": n, "PopTotal": pop_total,
            "p-value": pv,
            "Enrichment_score": (k / list_total) / (n / pop_total),
            "Substances": ",".join(sorted(genes)),
        })
    df = pd.DataFrame(rows).sort_values("p-value")
    p = df["p-value"].values
    order = np.argsort(p)
    adj = p[order] * len(p) / (np.arange(len(p)) + 1)
    adj = np.clip(np.minimum.accumulate(adj[::-1])[::-1], 0, 1)
    q = np.empty(len(p)); q[order] = adj
    df["q-value"] = q
    return df[["id", "Classification_level1", "Classification_level2", "Term",
               "ListHits", "ListTotal", "PopHits", "PopTotal",
               "p-value", "q-value", "Enrichment_score", "Substances"]]


def main():
    sd = Path(__file__).resolve().parent
    f_counts = find_file(sd, ["host_counts_integer_for_DESeq2.tsv", "host_counts.tsv"],
                         desc="宿主 counts 矩阵")
    f_anno = find_file(sd, ["counts_anno.xls", "counts_anno.txt"], desc="注释表")
    f_fa = find_file(sd, ["Unigene.fa", "Unigene.fasta"], required=False)

    print(f"宿主矩阵：{f_counts}")
    print(f"注释表  ：{f_anno}")
    print(f"序列文件：{f_fa if f_fa else '未找到，FPKM 将以 CPM 近似并加注'}")
    print()

    host = pd.read_csv(f_counts, sep="\t")
    idcol = host.columns[0]
    samples = [c for c in host.columns if re.fullmatch(r"(H|S1|S2)-\d+", str(c))]
    host = host.set_index(idcol)[samples]
    print(f"宿主 Unigene：{len(host)}　样本：{len(samples)}")

    keep = host.mean(axis=1) > MIN_MEAN_COUNT
    host = host[keep]
    print(f"counts 均值 > {MIN_MEAN_COUNT} 后保留：{len(host)}")

    meta = pd.DataFrame({
        "stage": [s.split("-")[0] for s in samples],
        "tree": ["T" + s.split("-")[1] for s in samples],
    }, index=samples)
    print(f"配对模型：design = ~ tree + stage（单株 {meta['tree'].nunique()} 棵）")
    print()

    print("正在运行 DESeq2（pydeseq2），约需 1–3 分钟…")
    dds = DeseqDataSet(counts=host.T.astype(int), metadata=meta,
                       design="~tree + stage", refit_cooks=True, quiet=True)
    dds.deseq2()
    normed = pd.DataFrame(dds.layers["normed_counts"], index=samples,
                          columns=host.index).T
    print("完成。\n")

    # FPKM：宿主子集内重新归一化
    if f_fa:
        lens = read_fasta_lengths(f_fa)
        kb = host.index.map(lambda x: lens.get(x, np.nan)) / 1000.0
        libs = host.sum(axis=0) / 1e6
        fpkm = host.div(pd.Series(kb, index=host.index), axis=0).div(libs, axis=1)
        fpkm_note = "宿主子集内重归一化 FPKM"
    else:
        libs = host.sum(axis=0) / 1e6
        fpkm = host.div(libs, axis=1)
        fpkm_note = "CPM（未找到序列文件，无法按长度归一）"

    anno = pd.read_csv(f_anno, sep="\t", low_memory=False)
    anno = anno.loc[:, [c for c in anno.columns if not str(c).startswith("Unnamed:")]]
    acols = [c for c in ANNO_COLS if c in anno.columns]
    anno = anno.set_index(anno.columns[0])[acols]

    host_set = set(host.index)
    pairs, pw_name, pw_cls, kegg_source = load_kegg(sd, f_anno, host_set)
    bg_pairs = [(g, k) for g, k in pairs if g in host_set]
    n_bg = len({g for g, _ in bg_pairs})
    n_pw = len({k for _, k in bg_pairs})
    print(f"KEGG 注释来源：{kegg_source}")
    print(f"宿主背景中带通路注释的 Unigene：{n_bg} / {len(host)} "
          f"({100*n_bg/len(host):.1f}%)，覆盖 {n_pw} 条通路")
    print(f"基因-通路对：{len(bg_pairs)}（一基因可对多通路）")
    print()

    out = sd / "2.gene"
    (out / "2.2.different_expressed_gene").mkdir(parents=True, exist_ok=True)
    (out / "2.3.function_of_different_expressed_gene" / "2.KEGG_enrichment").mkdir(
        parents=True, exist_ok=True)

    summary = []
    for case, ctrl in CONTRASTS:
        tag = f"{case}-vs-{ctrl}"
        print(f"— {tag} —")
        st = DeseqStats(dds, contrast=["stage", case, ctrl], quiet=True)
        st.summary()
        r = st.results_df.reindex(host.index)

        c_case = [s for s in samples if s.startswith(case + "-")]
        c_ctrl = [s for s in samples if s.startswith(ctrl + "-")]

        df = pd.DataFrame(index=host.index)
        df["gene_id"] = host.index
        df["BaseMean_DESeq2"] = r["baseMean"]
        df[f"BaseMean_control_{ctrl}"] = normed[c_ctrl].mean(axis=1)
        df[f"BaseMean_case_{case}"] = normed[c_case].mean(axis=1)
        df["FoldChange"] = 2 ** r["log2FoldChange"]
        df["log2FoldChange"] = r["log2FoldChange"]
        df["p-value"] = r["pvalue"]
        df["q-value"] = r["padj"]
        df["Regulation"] = np.where(r["log2FoldChange"] >= 0, "Up", "Down")
        for s in c_case + c_ctrl:
            df[f"FPKM_{s}"] = fpkm[s]
        for c in acols:
            df[c] = anno[c].reindex(host.index).fillna("--")

        d1 = out / "2.2.different_expressed_gene" / tag
        d1.mkdir(exist_ok=True)
        df.to_csv(d1 / f"{tag}-all.gene.xls", sep="\t", index=False)

        sig = df[(df["q-value"] < Q_CUT) &
                 (df["log2FoldChange"].abs() >= np.log2(FC_CUT))]
        sig.to_csv(d1 / f"{tag}-diff-q-val-{Q_CUT}-FC-{int(FC_CUT)}.gene.xls",
                   sep="\t", index=False)

        up = int((sig["Regulation"] == "Up").sum())
        dn = int((sig["Regulation"] == "Down").sum())
        print(f"  差异 Unigene：{len(sig)}（上调 {up} / 下调 {dn}）")

        e = enrich(set(sig["gene_id"]), pairs, host_set, pw_name, pw_cls)
        d2 = (out / "2.3.function_of_different_expressed_gene" /
              "2.KEGG_enrichment" / tag)
        d2.mkdir(exist_ok=True)
        if len(e):
            e.to_csv(d2 / f"enrichment-kegg-{tag}-Total.xls", sep="\t", index=False)
            ns = int((e["q-value"] < 0.05).sum())
            print(f"  KEGG 富集：{len(e)} 条通路，q < 0.05 者 {ns} 条")
            for _, row in e.head(5).iterrows():
                print(f"    {row['id']}  {str(row['Term'])[:44]:<44} "
                      f"{row['ListHits']:>3}/{row['ListTotal']:<4} q={row['q-value']:.3g}")
        else:
            ns = 0
            print("  KEGG 富集：无结果（差异基因中无带通路注释者）")
        summary.append((tag, len(sig), up, dn, len(e), ns))
        print()

    lines = ["=" * 62, "宿主子集差异分析 + 宿主背景 KEGG 富集", "=" * 62,
             f"宿主 Unigene（过滤后）: {len(host)}",
             f"模型: ~ tree + stage（配对）",
             f"差异阈值: q < {Q_CUT} 且 |FoldChange| >= {FC_CUT}",
             f"FPKM 口径: {fpkm_note}",
             f"KEGG 注释来源: {kegg_source}",
             f"KEGG 背景: {n_bg} 条宿主 Unigene，{n_pw} 条通路，{len(bg_pairs)} 个基因-通路对", "",
             f"{'比较':<10}{'差异数':>8}{'上调':>8}{'下调':>8}{'通路数':>8}{'q<0.05':>8}"]
    for t, n, u, d, np_, ns in summary:
        lines.append(f"{t:<10}{n:>8}{u:>8}{d:>8}{np_:>8}{ns:>8}")
    lines += ["", "输出目录：" + str(out),
              "  可直接将 2.gene 文件夹打包为 zip 上传欧易云联合分析", "",
              "⚠ 局限：针叶树特化的树脂酸合成酶（左旋海松二烯合酶、松脂二烯醇氧化酶、",
              "   CYP720B 家族）在 KEGG 参考通路中无直系同源记录，本研究 13 条相关宿主",
              "   Unigene 无一获得 KO 编号，故二萜树脂酸合成不可能出现在富集结果中。",
              "   该通路的结论须以通路内多基因方向一致性与终产物鉴定为依据。", "=" * 62]
    text = "\n".join(lines)
    print(text)
    (sd / "宿主差异与富集摘要.txt").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
        print("\n完成。")
    except Exception as e:
        print()
        print("=" * 62)
        print("运行出错：")
        print(f"  {type(e).__name__}: {e}")
        print()
        import traceback
        traceback.print_exc()
        print("=" * 62)
        print("把上面的报错内容整段复制给我，我来定位问题。")
    try:
        input("按回车键退出...")
    except EOFError:
        pass
