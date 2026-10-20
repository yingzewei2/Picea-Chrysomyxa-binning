@echo off
cd /d "%~dp0"
title Host DE + KEGG enrichment

echo ============================================================
echo   宿主子集差异分析（配对模型） + 宿主背景 KEGG 富集
echo ------------------------------------------------------------
echo   请把以下文件放在本 bat 所在的文件夹里：
echo       host_counts_integer_for_DESeq2.tsv
echo       counts_anno.xls
echo       Unigene.KEGG.gene.anno.xls
echo       Unigene.KEGG.classification.xls
echo       Unigene.fa  (可选)
echo.
echo   运行约需 5-10 分钟，中途会长时间无输出，属正常，请勿关闭。
echo ============================================================
echo.

where python >nul 2>nul
if errorlevel 1 (
    echo [错误] 找不到 python 命令，请确认已安装并勾选 Add Python to PATH
    pause
    exit /b 1
)

python "%~dp002_host_DE_KEGG.py"
pause
