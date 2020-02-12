#!/bin/bash
# 部署到 github pages 脚本
# 错误时终止脚本
set -e

# Commit changes.
rm -rf public

msg0="modifying site `date`"
if [ $# -eq 1 ]
  then msg0="$1"
fi

msg1="building site `date`"
if [ $# -eq 1 ]
  then msg1="$1"
fi

hugo -t even # if using a theme, replace with `hugo -t <YOURTHEME>`


git add -A
git commit -m "$msg0"

git push -f git@github.com:zuozqi/blog.git master

# 删除打包文件夹

# 打包。even 是主题

# 进入打包文件夹
cd public

# Add changes to git.

git init
git add -A
git commit -m "$msg1"

# 推送到githu
# nusr.github.io 只能使用 master分支
git push -f git@github.com:nusr/zuozqi.github.io.git master

# 回到原文件夹
cd ..