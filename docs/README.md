# PyAcoustiX 文档

这是PyAcoustiX/SAcouS项目的Sphinx文档。

## 文档结构

```
docs/
├── source/                 # 文档源文件
│   ├── index.rst          # 主文档入口
│   ├── user_guide/        # 用户指南
│   ├── api_reference/     # API参考
│   ├── theory/            # 理论基础
│   ├── examples/          # 示例
│   └── developer_guide/   # 开发者指南
├── build/                 # 构建输出
├── Makefile              # 构建脚本
└── requirements.txt      # 文档依赖
```

## 本地构建

### 安装依赖

```bash
pip install -r requirements.txt
```

### 构建HTML文档

```bash
cd docs
make html
```

构建后的文档将在 `docs/build/html/` 目录中。

### 查看文档

```bash
# 使用Python HTTP服务器
cd docs/build/html
python -m http.server 8000

# 然后在浏览器中打开 http://localhost:8000
```

## 在线文档

文档会自动部署到GitHub Pages:

- 地址: https://shaoqigit.github.io/PyXfem/

## 文档内容

### 用户指南
- 安装说明
- 快速开始
- 1D/2D/3D教程
- 命令行使用

### API参考
- Mesh模块
- Materials模块
- acxfem (FEM核心)
- acxtmm (传递矩阵法)
- acxmor (模态降阶)
- interface (接口)

### 理论基础
- 有限元方法基础
- Helmholtz方程
- JCA模型
- Biot理论
- 传递矩阵法

### 开发者指南
- 架构设计
- 贡献指南
- 测试指南

## 更新文档

1. 编辑 `docs/source/` 目录下的 `.rst` 文件
2. 运行 `make html` 构建
3. 提交更改到GitHub
4. GitHub Actions会自动部署到Pages

## 文档格式

使用reStructuredText格式:

- 标题使用下划线
- 代码块使用 `.. code-block::`
- 数学公式使用 `.. math::`
- 交叉引用使用 `:doc:` 和 `:ref:`

## 帮助

如有问题，请联系:
- 邮箱: shaoqiwu@outlook.com
- GitHub Issues: https://github.com/Shaoqigit/PyXfem/issues
