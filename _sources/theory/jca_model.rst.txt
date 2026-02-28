JCA模型
=======

Johnson-Champoux-Allard (JCA) 模型是描述刚性骨架多孔吸声材料的标准模型。

多孔材料简介
------------

多孔吸声材料具有复杂的微观结构:

- 固体骨架 (不可动)
- 连通的孔隙网络
- 孔隙中充满流体 (通常是空气)

声波在多孔材料中传播时，由于粘滞和热损耗，声能被转化为热能。

等效流体假设
------------

JCA模型采用**等效流体**假设:

- 固体骨架完全刚性
- 只考虑孔隙中流体的运动
- 用等效复密度和复体积模量描述宏观行为

等效密度
--------

考虑粘滞效应的Johnson等人模型:

.. math::

   \tilde{\rho}_{eq} = \frac{\rho_0 \alpha_\infty}{\phi} \left[ 1 + \frac{\sigma \phi}{j \omega \rho_0 \alpha_\infty} \sqrt{1 + \frac{4 j \omega \rho_0 \alpha_\infty^2 \eta}{\sigma^2 \phi^2 \Lambda^2}} \right]

其中:

- :math:`\rho_0`: 流体密度 [kg/m³]
- :math:`\phi`: 孔隙率 [-]
- :math:`\sigma`: 流阻 [Pa·s/m²]
- :math:`\alpha_\infty`: 高频曲折度 [-]
- :math:`\eta`: 动力粘度 [Pa·s]
- :math:`\Lambda`: 粘滞特征长度 [m]

等效体积模量
------------

考虑热效应的Champoux-Allard模型:

.. math::

   \tilde{K}_{eq} = \frac{\gamma P_0}{\phi} \left[ \gamma - \frac{\gamma - 1}{1 + \frac{8 \eta}{j \omega \rho_0 \Pr \Lambda'^2} \sqrt{1 + \frac{j \omega \rho_0 \Pr \Lambda'^2}{16 \eta}}} \right]^{-1}

其中:

- :math:`\gamma`: 比热比 (空气约1.4)
- :math:`P_0`: 大气压 [Pa]
- :math:`\Pr`: Prandtl数 (空气约0.71)
- :math:`\Lambda'`: 热特征长度 [m]

等效声速和阻抗
--------------

从等效密度和体积模量，可以计算:

**等效声速**:

.. math::

   \tilde{c}_{eq} = \sqrt{\frac{\tilde{K}_{eq}}{\tilde{\rho}_{eq}}}

**等效特性阻抗**:

.. math::

   \tilde{Z}_{eq} = \sqrt{\tilde{\rho}_{eq} \tilde{K}_{eq}}

**复波数**:

.. math::

   \tilde{k} = \omega \sqrt{\frac{\tilde{\rho}_{eq}}{\tilde{K}_{eq}}}

物理意义
--------

**孔隙率** (:math:`\phi`)

孔隙体积与总体积之比:

.. math::

   \phi = \frac{V_{pores}}{V_{total}}

典型值: 0.9 - 0.99 (泡沫材料)

**流阻** (:math:`\sigma`)

描述流体流经材料时的阻力:

.. math::

   \Delta P = \sigma \cdot v \cdot d

其中 :math:`\Delta P` 是压降，:math:`v` 是流速，:math:`d` 是厚度。

典型值: 10³ - 10⁶ Pa·s/m²

**曲折度** (:math:`\alpha_\infty`)

描述孔隙路径的复杂程度:

.. math::

   \alpha_\infty = \left( \frac{\text{实际路径长度}}{\text{直线路径长度}} \right)^2

典型值: 1.0 - 3.0

**特征长度**

- **粘滞特征长度** :math:`\Lambda`: 与粘滞损耗相关的孔隙尺寸
- **热特征长度** :math:`\Lambda'`: 与热损耗相关的孔隙尺寸

对于圆柱形孔隙: :math:`\Lambda = \Lambda' = R` (半径)

对于一般多孔材料: :math:`\Lambda' \approx 2\Lambda`

频率特性
--------

**低频极限** (:math:`\omega \to 0`)

- 等效密度趋于无穷 (流体被"锁定")
- 等效体积模量趋于 :math:`\gamma P_0 / \phi`

**高频极限** (:math:`\omega \to \infty`)

- 等效密度: :math:`\tilde{\rho}_{eq} \to \rho_0 \alpha_\infty / \phi`
- 等效体积模量: :math:`\tilde{K}_{eq} \to \gamma P_0 / \phi`

**过渡频率**

- **粘滞特征频率**: :math:`\omega_v = \frac{\sigma \phi}{\rho_0 \alpha_\infty}`
- **热特征频率**: :math:`\omega_t = \frac{16 \nu}{\Lambda'^2}`

吸声系数
--------

对于厚度为 :math:`d` 的多孔材料层， backed by rigid wall:

**表面阻抗**:

.. math::

   Z_s = -j \tilde{Z}_{eq} \cot(\tilde{k} d)

**吸声系数**:

.. math::

   \alpha = 1 - \left| \frac{Z_s - Z_0}{Z_s + Z_0} \right|^2

其中 :math:`Z_0 = \rho_0 c_0` 是空气特性阻抗。

参数测量
--------

JCA参数可以通过实验测量:

1. **孔隙率**: 通过密度测量或压汞法
2. **流阻**: 通过气流实验 (ISO 9053)
3. **曲折度**: 通过超声传播时间测量
4. **特征长度**: 通过吸声系数反演或微观成像

PyAcoustiX实现
--------------

.. code-block:: python

   from SAcouS.Materials import EquivalentFluid
   
   # 定义JCA材料
   porous = EquivalentFluid(
       name='foam',
       phi=0.98,           # 孔隙率
       sigma=3750,         # 流阻
       alpha=1.17,         # 曲折度
       Lambda_prime=742e-6,  # 热特征长度
       Lambda=110e-6       # 粘滞特征长度
   )
   
   # 设置频率
   omega = 2 * np.pi * 1000
   porous.set_frequency(omega)
   
   # 获取等效属性
   print(f"等效密度: {porous.rho_eq_til}")
   print(f"等效体积模量: {porous.K_eq_til}")
   print(f"等效声速: {porous.c_eq_til}")

进一步阅读
----------

- Johnson, D.L. et al. "Theory of dynamic permeability and tortuosity in fluid-saturated porous media"
- Champoux, Y. and Allard, J.F. "Dynamic tortuosity and bulk modulus in air-saturated porous media"
- Allard, J.F. and Atalla, N. "Propagation of Sound in Porous Media"
