传递矩阵法 (TMM)
=================

传递矩阵法(Transfer Matrix Method, TMM)是分析多层声学结构的高效方法。

基本原理
--------

TMM基于以下假设:

1. **平面波**: 声波以平面波形式传播
2. **垂直入射**: 波垂直于层界面传播
3. **无耦合**: 各层之间只有压力和速度的连续性

单层传递矩阵
------------

考虑厚度为 :math:`d` 的单层材料，声波沿 :math:`x` 方向传播。

在层内，声压和质点速度可以表示为:

.. math::

   p(x) = A e^{-jkx} + B e^{+jkx}

.. math::

   v(x) = \frac{1}{Z} (A e^{-jkx} - B e^{+jkx})

其中:

- :math:`A`: 正向传播波幅值
- :math:`B`: 反向传播波幅值
- :math:`k`: 波数
- :math:`Z`: 特性阻抗

**状态向量**

定义状态向量:

.. math::

   \mathbf{S} = \begin{bmatrix} p \\ v \end{bmatrix}

**传递矩阵**

层两端的声场关系:

.. math::

   \mathbf{S}(d) = T \cdot \mathbf{S}(0)

其中传递矩阵:

.. math::

   T = \begin{bmatrix} \cos(kd) & -jZ\sin(kd) \\ -j\frac{\sin(kd)}{Z} & \cos(kd) \end{bmatrix}

多层结构
--------

对于 :math:`N` 层结构，总传递矩阵为各层矩阵的乘积:

.. math::

   T_{total} = T_N \cdot T_{N-1} \cdots T_2 \cdot T_1

注意矩阵乘法的顺序: 从最后一层到第一层。

边界条件
--------

**入射端** (:math:`x=0`)

.. math::

   \mathbf{S}(0) = \begin{bmatrix} p_i + p_r \\ v_i + v_r \end{bmatrix}

其中 :math:`p_i, v_i` 是入射波，:math:`p_r, v_r` 是反射波。

**透射端** (:math:`x=L`)

.. math::

   \mathbf{S}(L) = \begin{bmatrix} p_t \\ v_t \end{bmatrix}

对于半无限透射介质:

.. math::

   p_t = Z_t v_t

声学性能指标
------------

**反射系数**

.. math::

   R = \frac{p_r}{p_i} = \frac{Z_s - Z_0}{Z_s + Z_0}

其中:

- :math:`Z_0`: 入射介质阻抗
- :math:`Z_s`: 表面阻抗

**吸声系数**

.. math::

   \alpha = 1 - |R|^2 = \frac{4 \text{Re}(Z_s)/Z_0}{|1 + Z_s/Z_0|^2}

**透射系数**

.. math::

   T = \frac{p_t}{p_i}

**传声损失**

.. math::

   TL = 20 \log_{10} \left| \frac{1}{T} \right| \quad \text{[dB]}

表面阻抗
--------

从传递矩阵可以得到表面阻抗:

.. math::

   Z_s = \frac{T_{11} Z_t + T_{12}}{T_{21} Z_t + T_{22}}

其中 :math:`Z_t` 是透射端阻抗。

对于刚性背衬 (:math:`v=0`):

.. math::

   Z_s = \frac{T_{12}}{T_{22}}

材料模型
--------

**流体层**

.. math::

   k = \frac{\omega}{c}, \quad Z = \rho c

**多孔材料 (JCA模型)**

.. math::

   k = \omega \sqrt{\frac{\tilde{\rho}_{eq}}{\tilde{K}_{eq}}}, \quad Z = \sqrt{\tilde{\rho}_{eq} \tilde{K}_{eq}}

**弹性层**

对于薄板，需要考虑弯曲波:

.. math::

   Z = j\omega m \left( 1 - \frac{f^2}{f_c^2} \right)

其中 :math:`m` 是面密度，:math:`f_c` 是临界频率。

斜入射
------

对于入射角 :math:`\theta`，修正传递矩阵:

.. math::

   k_x = k \cos\theta

.. math::

   Z_\theta = \frac{Z}{\cos\theta}

吸声系数随入射角变化:

- 0° (垂直入射): 最大吸声
- 90° (掠入射): 吸声最小

扩散场吸声系数
~~~~~~~~~~~~~~

对于随机入射，平均吸声系数:

.. math::

   \alpha_{diff} = 2 \int_0^{\pi/2} \alpha(\theta) \sin\theta \cos\theta \, d\theta

频率特性
--------

**四分之一波长共振**

对于厚度为 :math:`d` 的单层材料，当:

.. math::

   d = \frac{\lambda}{4} = \frac{c}{4f}

时出现最大吸声。

**低频行为**

- 薄层: 吸声很小
- 增加厚度可以提高低频吸声

**高频行为**

- 吸声趋于稳定值
- 可能出现干涉效应

TMM vs FEM
----------

**TMM优势**:

- 计算极快 (矩阵乘法)
- 适合参数化研究
- 适合优化设计

**TMM局限**:

- 只适用于平面波
- 不考虑边缘效应
- 复杂几何困难

**FEM-TMM耦合**:

PyAcoustiX支持FEM-TMM耦合，结合两者的优势:

- 复杂几何区域用FEM
- 多层均匀区域用TMM

PyAcoustiX实现
--------------

.. code-block:: python

   from SAcouS.acxtmm import TMM
   from SAcouS.Materials import Air, EquivalentFluid
   
   # 定义多层结构
   air = Air('air')
   foam = EquivalentFluid('foam', 0.98, 3750, 1.17, 742e-6, 110e-6)
   
   materials = [air, foam, air]
   thicknesses = [0.01, 0.05, 0.01]  # m
   
   # 创建TMM模型
   freq = 1000
   tmm = TMM(materials, thicknesses, freq)
   
   # 计算性能指标
   alpha = tmm.absorption_coefficient()
   Zs = tmm.surface_impedance()
   TL = tmm.transmission_loss()
   
   print(f"吸声系数: {alpha:.3f}")
   print(f"表面阻抗: {Zs:.2f}")
   print(f"传声损失: {TL:.2f} dB")

进一步阅读
----------

- Beranek, L.L. and Ver, I.L. "Noise and Vibration Control Engineering"
- Cox, T.J. and D'Antonio, P. "Acoustic Absorbers and Diffusers"
- Allard, J.F. and Atalla, N. "Propagation of Sound in Porous Media"
