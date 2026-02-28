Materials模块
=============

Materials模块定义了各种声学材料模型，包括流体、多孔材料和弹性材料。

.. module:: SAcouS.Materials

基类
----

.. autoclass:: SAcouS.Materials.BaseMaterial
   :members:
   :undoc-members:
   :show-inheritance:

流体材料
--------

空气
~~~~

.. autoclass:: SAcouS.Materials.Air
   :members:
   :undoc-members:
   :show-inheritance:

通用流体
~~~~~~~~

.. autoclass:: SAcouS.Materials.Fluid
   :members:
   :undoc-members:
   :show-inheritance:

多孔材料
--------

JCA等效流体
~~~~~~~~~~~

.. autoclass:: SAcouS.Materials.EquivalentFluid
   :members:
   :undoc-members:
   :show-inheritance:

Limp多孔材料
~~~~~~~~~~~~

.. autoclass:: SAcouS.Materials.LimpPorousMaterial
   :members:
   :undoc-members:
   :show-inheritance:

弹性材料
--------

.. autoclass:: SAcouS.Materials.ElasticMaterial
   :members:
   :undoc-members:
   :show-inheritance:

Biot多孔弹性材料
----------------

.. autoclass:: SAcouS.Materials.PoroElasticMaterial
   :members:
   :undoc-members:
   :show-inheritance:

材料工厂
--------

.. autoclass:: SAcouS.Materials.MaterialFactory
   :members:
   :undoc-members:
   :show-inheritance:

使用示例
--------

预定义材料
~~~~~~~~~~

.. code-block:: python

   from SAcouS.Materials import Air
   
   # 标准空气
   air = Air('air')
   print(f"密度: {air.rho} kg/m³")
   print(f"声速: {air.c} m/s")

自定义流体
~~~~~~~~~~

.. code-block:: python

   from SAcouS.Materials import Fluid
   
   # 自定义流体 (水)
   water = Fluid('water', rho=1000, c=1500)

JCA多孔材料
~~~~~~~~~~~

.. code-block:: python

   from SAcouS.Materials import EquivalentFluid
   
   # 定义JCA参数
   phi = 0.98           # 孔隙率
   sigma = 3750         # 流阻 [Pa·s/m²]
   alpha = 1.17         # 曲折度
   Lambda_prime = 742e-6  # 热特征长度 [m]
   Lambda = 110e-6      # 粘滞特征长度 [m]
   
   porous = EquivalentFluid('porous', phi, sigma, alpha, 
                            Lambda_prime, Lambda)
   
   # 设置频率 (计算等效属性)
   omega = 2 * np.pi * 1000
   porous.set_frequency(omega)
   
   print(f"等效密度: {porous.rho_eq_til}")
   print(f"等效体积模量: {porous.K_eq_til}")

Biot材料
~~~~~~~~

.. code-block:: python

   from SAcouS.Materials import PoroElasticMaterial
   
   # JCA参数 + 弹性参数
   E = 1e6          # 杨氏模量 [Pa]
   nu = 0.3         # 泊松比
   eta = 0.01       # 损耗因子
   rho_solid = 30   # 固体密度 [kg/m³]
   
   biot = PoroElasticMaterial('biot', phi, sigma, alpha,
                               Lambda_prime, Lambda, rho_solid,
                               E, nu, eta)

材料工厂
~~~~~~~~

.. code-block:: python

   from SAcouS.Materials import MaterialFactory
   
   # 使用工厂创建材料
   air = MaterialFactory.create_material('AIR', 'air')
   porous = MaterialFactory.create_material('RIGID_POROUS', 'foam',
                                            0.98, 3750, 1.17, 
                                            742e-6, 110e-6)
