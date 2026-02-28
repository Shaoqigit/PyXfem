Biot理论
========

Biot理论描述了孔隙弹性材料(可变形骨架)中的波传播，是JCA模型的扩展。

Biot理论概述
------------

Maurice Biot在1956年提出了描述饱和孔隙介质中波传播的理论。与JCA模型不同，Biot理论:

- 考虑固体骨架的弹性变形
- 考虑固体与流体之间的惯性耦合
- 预测两种压缩波和一种剪切波

基本假设
--------

1. **小变形**: 线弹性理论适用
2. **连通孔隙**: 流体可以自由流动
3. **可压缩流体**: 流体具有可压缩性
4. **统计各向同性**: 材料宏观上各向同性
5. **可忽略热效应**: 等温过程

双相系统
--------

Biot理论将孔隙弹性材料视为双相系统:

- **固相** (骨架): 体积分数 :math:`1-\phi`
- **液相** (孔隙流体): 体积分数 :math:`\phi`

定义平均位移:

- :math:`\mathbf{u}`: 固体骨架位移
- :math:`\mathbf{U}`: 流体位移

相对位移:

.. math::

   \mathbf{w} = \phi (\mathbf{U} - \mathbf{u})

表示流体相对于固体的运动。

本构关系
--------

应力-应变关系:

.. math::

   \sigma_{ij} = (P - 2N) e \delta_{ij} + 2N e_{ij} + Q \epsilon \delta_{ij}

.. math::

   s = Q e + R \epsilon

其中:

- :math:`\sigma_{ij}`: 总应力张量
- :math:`s`: 流体压力
- :math:`e_{ij} = \frac{1}{2}(\partial_i u_j + \partial_j u_i)`: 固体应变
- :math:`e = \nabla \cdot \mathbf{u}`: 固体体积应变
- :math:`\epsilon = \nabla \cdot \mathbf{U}`: 流体体积应变

弹性系数
--------

Biot弹性系数 :math:`P, Q, R, N` 可以用更易测量的参数表示:

**剪切模量**:

.. math::

   N = \frac{E_s}{2(1+\nu_s)}

其中 :math:`E_s` 和 :math:`\nu_s` 是固体骨架的杨氏模量和泊松比。

**其他系数**:

.. math::

   P = \frac{4}{3}N + K_b + \frac{(1-\phi)^2}{\phi}K_f

.. math::

   Q = (1-\phi)K_f

.. math::

   R = \phi K_f

其中:

- :math:`K_b`: 骨架体积模量
- :math:`K_f`: 流体体积模量
- :math:`K_s`: 固体颗粒体积模量

运动方程
--------

**固体运动方程**:

.. math::

   \rho_{11} \ddot{u}_i + \rho_{12} \ddot{U}_i + b(\dot{u}_i - \dot{U}_i) = \frac{\partial \sigma_{ij}}{\partial x_j}

**流体运动方程**:

.. math::

   \rho_{12} \ddot{u}_i + \rho_{22} \ddot{U}_i - b(\dot{u}_i - \dot{U}_i) = \frac{\partial s}{\partial x_i}

其中:

- :math:`\rho_{11} = (1-\phi)\rho_s - \rho_{12}`: 固体有效密度
- :math:`\rho_{22} = \phi\rho_f - \rho_{12}`: 流体有效密度
- :math:`\rho_{12}`: 惯性耦合系数
- :math:`b = \frac{\sigma \phi^2}{\omega}`: 粘滞耦合系数

波类型
------

Biot理论预测三种体波:

**快压缩波 (P1)**

- 固体和流体同相运动
- 速度接近纯固体中的波速
- 衰减较小

**慢压缩波 (P2)**

- 固体和流体反相运动
- 速度较慢
- 衰减很大，传播距离短

**剪切波 (S)**

- 横波，只有固体骨架参与
- 流体不影响剪切波速

波速公式
--------

对于简谐波 :math:`e^{j(\omega t - kx)}`，色散关系给出:

.. math::

   \delta_{1,2}^2 = \frac{\omega^2}{2} \left[ \frac{\delta_s^2 + \delta_f^2 \pm \sqrt{(\delta_s^2 - \delta_f^2)^2 + 4\gamma^2 \delta_f^2 \delta_s^2}}{1 - \gamma^2} \right]

其中:

- :math:`\delta_1`: 快波波数
- :math:`\delta_2`: 慢波波数
- :math:`\delta_s = \omega\sqrt{\rho_s/P}`: 纯固体波数
- :math:`\delta_f = \omega\sqrt{\rho_f/K_f}`: 纯流体波数
- :math:`\gamma = \frac{Q}{R}\sqrt{\frac{R}{P}}`: 耦合系数

边界条件
--------

**开放边界** (流体自由流动):

.. math::

   s = 0

**封闭边界** (不可渗透):

.. math::

   w_n = 0

**滑动边界** (无摩擦):

.. math::

   \tau_{nt} = 0

其中 :math:`\tau_{nt}` 是切应力。

频率特性
--------

**低频极限** (:math:`\omega \to 0`)

- 流体被锁定，表现为单相材料
- 有效密度: :math:`\rho_{eq} = (1-\phi)\rho_s + \phi\rho_f`

**高频极限** (:math:`\omega \to \infty`)

- 惯性耦合主导
- 快波和慢波解耦

**特征频率**

- **粘滞特征频率**: :math:`f_v = \frac{\sigma \phi}{2\pi \rho_f}`
- 在此频率附近，慢波开始出现

应用
----

Biot理论应用于:

- **地球物理**: 地震波在饱和岩石中的传播
- **生物力学**: 骨组织、软组织的声学特性
- **环境声学**: 饱和土壤中的噪声传播
- **石油工程**: 油藏声学测井

PyAcoustiX实现
--------------

.. code-block:: python

   from SAcouS.Materials import PoroElasticMaterial
   
   # 定义Biot材料
   biot = PoroElasticMaterial(
       name='saturated_soil',
       phi=0.4,            # 孔隙率
       sigma=1e4,          # 流阻
       alpha=2.0,          # 曲折度
       Lambda_prime=100e-6,  # 热特征长度
       Lambda=50e-6,       # 粘滞特征长度
       rho_solid=2650,     # 固体密度
       E=5e9,              # 杨氏模量
       nu=0.3,             # 泊松比
       eta=0.01            # 损耗因子
   )
   
   # 设置频率
   omega = 2 * np.pi * 500
   biot.set_frequency(omega)
   
   # 获取波数
   print(f"快波波数: {biot.delta_1}")
   print(f"慢波波数: {biot.delta_2}")
   print(f"剪切波波数: {biot.delta_3}")

进一步阅读
----------

- Biot, M.A. "Theory of propagation of elastic waves in a fluid-saturated porous solid"
- Bourbie, T. et al. "Propagation of elastic waves in porous media"
- Carcione, J.M. "Wave Fields in Real Media"
