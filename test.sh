#!/bin/bash
echo "=== Alice O2编译环境诊断 ==="
echo "时间: $(date)"
echo

echo "1. 系统信息:"
sw_vers
sysctl -n machdep.cpu.brand_string
echo "内存: $(sysctl -n hw.memsize | awk '{print $1/1024/1024/1024" GB"}')"
echo

echo "2. Xcode状态:"
xcode-select -p
xcodebuild -version 2>/dev/null || echo "Xcode未安装"
ls -la /Applications/Xcode.app 2>/dev/null || echo "Xcode.app不存在"
echo

echo "3. 编译器状态:"
for cmd in clang clang++ gcc g++; do
    which $cmd 2>/dev/null && $cmd --version | head -1
done
echo

echo "4. aliBuild状态:"
aliBuild version 2>/dev/null || echo "aliBuild未找到"
echo "ALIBUILD_WORK_DIR: $ALIBUILD_WORK_DIR"
echo

echo "5. 磁盘状态:"
df -h /Users/allenzhou/ALICE/sw
echo

echo "6. 进程限制:"
ulimit -a | grep -E "(open files|max user processes)"
echo

echo "7. 失败的构建目录:"
ls -la /Users/allenzhou/ALICE/sw/BUILD/ 2>/dev/null | wc -l | xargs echo "构建目录数:"
echo

echo "8. CMake状态:"
cmake --version 2>/dev/null | head -1 || echo "CMake未找到"
