@echo off
cd /d "%~dp0"
title 转录组宿主/病原分箱

echo ============================================================
echo   转录组宿主 / 病原 分箱
echo ------------------------------------------------------------
echo   请把 counts_anno.xls 放在本 bat 所在的文件夹里，
echo   然后直接运行本文件即可。
echo ============================================================
echo.

where python >nul 2>nul
if errorlevel 1 (
    echo [错误] 系统里找不到 python 命令。
    echo.
    echo 可能的原因：
    echo   1. 没有安装 Python  -^> 去 https://www.python.org/downloads/ 下载
    echo   2. 安装时没有勾选 "Add Python to PATH"  -^> 重新安装并勾选该项
    echo.
    pause
    exit /b 1
)

python "%~dp001_binning_host_pathogen.py"

if errorlevel 1 (
    echo.
    echo [提示] 脚本以非零状态退出，请把上面的报错内容整段复制出来。
    pause
)
