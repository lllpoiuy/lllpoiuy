# 数学分析（二）后半学期-复习笔记

本文为 deepseek 辅助复习，希望能够帮助到希望复习 de Rham 上同调群、同伦群、同调群、Gauss-Bonnet 等的同学。



## 0. 微分流形 (Differentiable Manifolds)

* **定义：** 一个 **n 维微分流形** $M$ 是一个满足以下条件的拓扑空间：

  1.  **局部欧氏性 (Locally Euclidean):** 每一点 $p \in M$ 都有一个邻域 $U$ 同胚于 $\mathbb{R}^n$ 的一个开子集。
  2.  **Hausdorff 与第二可数性 (Hausdorff and Second Countable):** 分离公理确保点可分，第二可数性保证存在可数拓扑基。
  3.  **光滑结构 (Smooth Structure):** 存在一个 **光滑图册 (smooth atlas)** $\mathcal{A} = \{(U_\alpha, \phi_\alpha)\}$，其中：
      *   $U_\alpha$ 是覆盖 $M$ 的开集。
      *   $\phi_\alpha: U_\alpha \to \phi_\alpha(U_\alpha) \subset \mathbb{R}^n$ 是同胚映射（坐标卡）。
      *   对于任意两个坐标卡 $(U_\alpha, \phi_\alpha)$ 和 $(U_\beta, \phi_\beta)$，如果 $U_\alpha \cap U_\beta \neq \emptyset$，则 **坐标变换 (transition map)** $\phi_\beta \circ \phi_\alpha^{-1}: \phi_\alpha(U_\alpha \cap U_\beta) \to \phi_\beta(U_\alpha \cap U_\beta)$ 是光滑的 ($C^\infty$) 映射。

  *   这个光滑图册定义了流形 $M$ 上 **光滑函数 (smooth functions)** 的概念：一个函数 $f: M \to \mathbb{R}$ 是光滑的，如果对于每个坐标卡 $(U_\alpha, \phi_\alpha)$，函数 $f \circ \phi_\alpha^{-1}: \phi_\alpha(U_\alpha) \to \mathbb{R}$ 是光滑的。



## 1. 微分形式 (Differential Forms)

### 1.1 定义在局部坐标卡上的微分形式

在流形 $M$ 上，我们可以定义更复杂的几何对象：**微分形式**。它们是定义在切丛和余切丛上的张量场的反对称化。

* **余切丛 (Cotangent Bundle):** 在每一点 $p \in M$，切空间 $T_pM$ 的对偶空间称为余切空间 $T_p^*M$。$T_p^*M$ 中的元素称为 **余向量 (covectors)** 或 **1-形式 (1-forms)** 在点 $p$。所有余切空间的并集 $\bigcup_{p\in M} T_p^*M$ 构成余切丛 $T^*M$。

* **k-形式 (k-forms):** 一个 **k 次微分形式 (differential k-form)** $\omega$ 是一个光滑的、反对称的协变 k-张量场。具体来说，它在每一点 $p$ 指派了一个反对称的多线性映射：
  $$
  \omega_p: \underbrace{T_pM \times \cdots \times T_pM}_{k\ \text{times}} \to \mathbb{R}
  $$
  并且这个指派在坐标变换下光滑变化。

* **外代数 (Exterior Algebra):** 在一点 $p$ 的所有 k-形式构成一个向量空间 $\bigwedge^k T_p^*M$。所有点上的 k-形式的光滑截面构成一个向量空间，记为 $\Omega^k(M)$。

  *   $\Omega^0(M) = C^\infty(M)$：光滑函数（0-形式）。
  *   $\Omega^1(M)$：1-形式（余向量场）。
  *   $\Omega^k(M) = 0$ 当 $k > \dim M$ 或 $k < 0$。

* **外积 (Wedge Product):** 我们可以定义一种乘法运算 $\wedge$（外积）将形式乘在一起：
  $$
  \wedge: \Omega^k(M) \times \Omega^l(M) \to \Omega^{k+l}(M)
  $$
  这个运算满足 **反对称性 (antisymmetry)**：$\alpha \wedge \beta = (-1)^{kl} \beta \wedge \alpha$（$\alpha$ 是 k-形式，$\beta$ 是 l-形式），并且是双线性和结合律的（在结合律的意义上）。

* **外微分 (Exterior Derivative):** 这是 de Rham 上同调的核心算子。它是一个线性算子：
  $$
  d: \Omega^k(M) \to \Omega^{k+1}(M)
  $$
  满足以下关键性质：

  1. **d 在函数上的作用：** 如果 $f \in \Omega^0(M) = C^\infty(M)$，则 $df$ 是其微分，这是一个 1-形式。在局部坐标 $(x^1, ..., x^n)$ 下，$df = \sum_{i=1}^n \frac{\partial f}{\partial x^i} dx^i$。

  2. **d 在形式上的作用：** 如果 $\omega = \sum_I a_I dx^I \in \Omega^k(M)$（其中 $I = (i_1 < ... < i_k)$ 是多重指标，$dx^I = dx^{i_1} \wedge \cdots \wedge dx^{i_k}$)，则：
     $$
     d\omega = \sum_I da_I \wedge dx^I = \sum_I \sum_j \frac{\partial a_I}{\partial x^j} dx^j \wedge dx^I
     $$

  3. **莱布尼茨律 (Product Rule):** $d(\alpha \wedge \beta) = d\alpha \wedge \beta + (-1)^{\deg \alpha} \alpha \wedge d\beta$。

  4. **幂零性 (Nilpotence):** 最重要的性质！对任意形式 $\omega$，有：
     $$
     d(d\omega) = d^2\omega = 0
     $$
     这个性质 $d^2 = 0$ 是整个上同调理论的基石。



### 1.2 定义在全局上的微分流形

+ 在微分流形上，微分形式在**局部坐标卡**上定义，但通过以下机制，它们能粘合成**全局**的微分形式空间 $\Omega^k(M)$：

  #### 1.2.1. **局部定义与坐标变换的相容性**

  设 $M$ 是 $n$ 维微分流形，$\{(U_\alpha, \phi_\alpha)\}$ 是其坐标图册。  

  - 在坐标卡 $(U_\alpha, \phi_\alpha)$ 上，$\phi_\alpha: U_\alpha \to \mathbb{R}^n$ 给出局部坐标 $(x^1, \dots, x^n)$。  

  - 一个 $k$-形式 $\omega$ 在 $U_\alpha$ 上可表示为：
    $$
    \omega|_{U_\alpha} = \sum_{1 \leq i_1 < \cdots < i_k \leq n} f_{i_1 \cdots i_k}^\alpha  dx^{i_1} \wedge \cdots \wedge dx^{i_k}
    $$
    其中 $f_{i_1 \cdots i_k}^\alpha: U_\alpha \to \mathbb{R}$ 是光滑函数。

  **关键问题**：当两个坐标卡 $U_\alpha$ 和 $U_\beta$ 相交（即 $U_\alpha \cap U_\beta \neq \emptyset$）时，如何保证 $\omega$ 在交集上一致？

  #### 1.2.2. **微分形式的坐标变换规则**

  设 $U_\alpha \cap U_\beta \neq \emptyset$，坐标变换 $\phi_\beta \circ \phi_\alpha^{-1}: \phi_\alpha(U_\alpha \cap U_\beta) \to \phi_\beta(U_\alpha \cap U_\beta)$ 是光滑映射。  

  - 若 $\phi_\alpha$ 给出坐标 $x$，$\phi_\beta$ 给出坐标 $y$，则 $y = y(x)$。  

  - 微分形式的变换由 **Jacobian 矩阵** 控制：
    $$
    dy^j = \sum_{i=1}^n \frac{\partial y^j}{\partial x^i}  dx^i
    $$

  - 对 $k$-形式，其变换规则为：
    $$
    dy^{j_1} \wedge \cdots \wedge dy^{j_k} = \sum_{1 \leq i_1 < \cdots < i_k \leq n} \det\left( \frac{\partial y^{j_\bullet}}{\partial x^{i_\bullet}} \right)  dx^{i_1} \wedge \cdots \wedge dx^{i_k}
    $$
    其中 $\det\left( \frac{\partial y^{j_\bullet}}{\partial x^{i_\bullet}} \right)$ 是 Jacobian 矩阵的 $k \times k$ 子式行列式。

  #### 1.2.3. **全局定义的相容性条件**

  为使 $\omega$ 在 $U_\alpha \cap U_\beta$ 上一致，需满足：
  $$
  \omega|_{U_\alpha} = \omega|_{U_\beta} \quad \text{在 } U_\alpha \cap U_\beta \text{ 上}
  $$
  即分量函数满足：
  $$
  f_{j_1 \cdots j_k}^\beta = \sum_{1 \leq i_1 < \cdots < i_k \leq n} f_{i_1 \cdots i_k}^\alpha \cdot \det\left( \frac{\partial x^{i_\bullet}}{\partial y^{j_\bullet}} \right)
  $$
  此条件确保：

  - $\omega$ 的定义不依赖于坐标卡的选择。
  - 局部表达式能唯一粘合成整体流形上的微分形式。

  > **示例**：圆周 $S^1$ 上的 $d\theta$  
  >
  > - 在坐标卡 $U_1 = S^1 \setminus \{(1,0)\}$ 中，$\theta_1 \in (0, 2\pi)$，定义 $\omega|_{U_1} = d\theta_1$。  
  > - 在坐标卡 $U_2 = S^1 \setminus \{(-1,0)\}$ 中，$\theta_2 \in (-\pi, \pi)$，定义 $\omega|_{U_2} = d\theta_2$。  
  > - 在 $U_1 \cap U_2$ 上，$\theta_1 = \theta_2 + c$（常数偏移），故 $d\theta_1 = d\theta_2$。  
  >   因此 $\omega$ 是全局定义的 1-形式。

  #### 1.2.4. **为什么 $\Omega^k(M)$ 是全局空间？**

  $\Omega^k(M)$ 定义为：
  $$
  \Omega^k(M) = \left\{ \omega: \omega \text{ 是 } M \text{ 上光滑的 } k\text{-形式} \right\}
  $$
  其全局性由以下保证：

  1. **局部到整体的粘合**：  
     通过坐标变换的相容性，局部形式可唯一粘合成整体形式（类似向量丛的截面）。
  2. **光滑结构的协调**：  
     流形的微分结构要求坐标变换光滑，这确保了微分形式的分量函数光滑过渡。
  3. **外代数的函子性**：  
     微分形式是反对称协变张量场，其变换规则由切丛与余切丛的拉回映射自然导出。

  #### 1.2.5. **与向量丛的关系**

  $\Omega^k(M)$ 本质是 **外 $k$-形式丛 $\bigwedge^k T^*M$** 的光滑截面空间：

  - 每点 $p \in M$ 的纤维 $\bigwedge^k T_p^*M$ 是一个向量空间。
  - 微分形式 $\omega$ 是截面：$p \mapsto \omega_p \in \bigwedge^k T_p^*M$。
  - 坐标变换规则保证了截面的整体光滑性。



### 1.3 外微分的运算性质、例子及与向量场运算的关系

#### 1.3.1 外微分的基本运算性质

设 $ M $ 为光滑流形，$ \omega, \eta $ 为微分形式，$ f $ 为光滑函数。外微分算子 $ d: \Omega^k(M) \to \Omega^{k+1}(M) $ 满足以下性质：

1. **线性性**  
   $$
   d(a\omega + b\eta) = a \, d\omega + b \, d\eta \quad (a,b \in \mathbb{R})
   $$

2. **Leibniz 法则（反导子性质）**  
   $$
   d(\omega \wedge \eta) = d\omega \wedge \eta + (-1)^{\deg \omega} \omega \wedge d\eta
   $$

3. **幂零性**  
   $$
   d^2 = d \circ d = 0
   $$

4. **函数的外微分**  
   对 0-形式（函数）$ f $，$ df $ 是余切向量场：  
   $$
   df = \sum_{i=1}^n \frac{\partial f}{\partial x^i} dx^i
   $$

#### 1.3.2 外微分计算示例（以 $ \mathbb{R}^3 $ 为例）

##### 例 1：0-形式（函数）的外微分

设 $ f(x,y,z) = x^2 y + \sin z $，则：  
$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy + \frac{\partial f}{\partial z}dz = 2xy \, dx + x^2 \, dy + \cos z \, dz
$$

##### 例 2：1-形式的外微分

设 $ \omega = yz \, dx + xz \, dy + xy \, dz $，则：  
$$
d\omega = d(yz) \wedge dx + d(xz) \wedge dy + d(xy) \wedge dz
$$
展开计算：  
$$
\begin{align*}
d(yz) &= z \, dy + y \, dz \\
d(xz) &= z \, dx + x \, dz \\
d(xy) &= y \, dx + x \, dy \\
\implies d\omega &= (z \, dy + y \, dz) \wedge dx + (z \, dx + x \, dz) \wedge dy + (y \, dx + x \, dy) \wedge dz \\
&= -\cancel{z \, dx \wedge dy} - y \, dx \wedge dz + z \, dx \wedge dy - \cancel{x \, dy \wedge dz} + y \, dx \wedge dz + \cancel{x \, dy \wedge dz} \\
&= 0 \quad (\text{所有项抵消})
\end{align*}
$$

##### 例 3：2-形式的外微分

设 $ \eta = x \, dy \wedge dz + y \, dz \wedge dx + z \, dx \wedge dy $，则：  
$$
d\eta = d(x) \wedge dy \wedge dz + d(y) \wedge dz \wedge dx + d(z) \wedge dx \wedge dy
$$
计算得：  
$$
d\eta = dx \wedge dy \wedge dz + dy \wedge dz \wedge dx + dz \wedge dx \wedge dy = 3 \, dx \wedge dy \wedge dz
$$
（因 $ dy \wedge dz \wedge dx = dx \wedge dy \wedge dz $ 等）

#### 1.3.3 外微分与向量场运算的关系（$ \mathbb{R}^3 $ 中）

在 $ \mathbb{R}^3 $ 中，微分形式与向量场通过以下对应关系等价：  

| 微分形式          | 向量场表示            | 记号                                                         |
| ----------------- | --------------------- | ------------------------------------------------------------ |
| 0-形式 $ f $      | 标量场 $ f $          | $ f $                                                        |
| 1-形式 $ \omega $ | 向量场 $ \mathbf{F} $ | $ \omega = F_x dx + F_y dy + F_z dz $                        |
| 2-形式 $ \eta $   | 向量场 $ \mathbf{G} $ | $ \eta = G_x dy \wedge dz + G_y dz \wedge dx + G_z dx \wedge dy $ |
| 3-形式 $ \nu $    | 标量场 $ h $          | $ \nu = h \, dx \wedge dy \wedge dz $                        |

##### 1.3.3.1. 梯度（Gradient）

- **向量场运算**：$ \nabla f = \left( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z} \right) $

- **外微分实现**：  
  $$
  d(f) = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy + \frac{\partial f}{\partial z}dz \quad \text{(1-形式)}
  $$

- **对应关系**：$ \nabla f \leftrightarrow df $

##### 1.3.3.2. 旋度（Curl）

- **向量场运算**：  
  $$
  \nabla \times \mathbf{F} = \left( \frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z}, \frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x}, \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right)
  $$

- **外微分实现**：  
  对 1-形式 $ \omega = F_x dx + F_y dy + F_z dz $:  
  $$
  d\omega = \left( \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right) dx \wedge dy + \left( \frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z} \right) dy \wedge dz + \left( \frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x} \right) dz \wedge dx
  $$

- **对应关系**：  
  $$
  \nabla \times \mathbf{F} \leftrightarrow d\omega \quad \text{(2-形式)}
  $$

##### 1.3.3.3. 散度（Divergence）

- **向量场运算**：  
  $$
  \nabla \cdot \mathbf{G} = \frac{\partial G_x}{\partial x} + \frac{\partial G_y}{\partial y} + \frac{\partial G_z}{\partial z}
  $$

- **外微分实现**：  
  对 2-形式 $ \eta = G_x dy \wedge dz + G_y dz \wedge dx + G_z dx \wedge dy $:  
  $$
  d\eta = \left( \frac{\partial G_x}{\partial x} + \frac{\partial G_y}{\partial y} + \frac{\partial G_z}{\partial z} \right) dx \wedge dy \wedge dz
  $$

- **对应关系**： 
  $$
  \nabla \cdot \mathbf{G} \leftrightarrow d\eta \quad \text{(3-形式)}
  $$

| 运算         | 向量场表示                       | 外微分操作                 | 阶数变化    |
| ------------ | -------------------------------- | -------------------------- | ----------- |
| **梯度**     | $ \nabla f $                     | $ d(f) $                   | $ 0 \to 1 $ |
| **旋度**     | $ \nabla \times \mathbf{F} $     | $ d(\omega_{\mathbf{F}}) $ | $ 1 \to 2 $ |
| **散度**     | $ \nabla \cdot \mathbf{G} $      | $ d(\eta_{\mathbf{G}}) $   | $ 2 \to 3 $ |
| **二阶零性** | $ \nabla \times (\nabla f) = 0 $ | $ d^2 f = 0 $              | -           |

外微分将梯度、旋度、散度统一为单一算子 $ d $，并在任意维流形上推广了向量分析，揭示了微分形式与拓扑的深刻联系。



## 2. de Rham 复形与上同调群

外微分算子 $d$ 和幂零性 $d^2=0$ 允许我们构造一个 **链复形 (chain complex)**，称为 **de Rham 复形 (de Rham complex)**：

$$
0 \longrightarrow \Omega^0(M) \xrightarrow{d} \Omega^1(M) \xrightarrow{d} \Omega^2(M) \xrightarrow{d} \cdots \xrightarrow{d} \Omega^n(M) \xrightarrow{d} 0 \xrightarrow{d} \cdots
$$

这个复形的核心特征是：一个映射 $d$ 的像（image）包含在下一个映射 $d$ 的核（kernel）中。这是因为 $d \circ d = 0$，意味着如果 $\beta = d\alpha$（$\beta$ 是恰当的），那么 $d\beta = d(d\alpha) = 0$（$\beta$ 也是闭的）。

* **闭形式 (Closed Forms):** 一个 k-形式 $\omega \in \Omega^k(M)$ 称为 **闭的 (closed)**，如果 $d\omega = 0$。所有闭 k-形式构成的向量空间是 $\ker(d: \Omega^k \to \Omega^{k+1})$，记为 $Z^k_{\text{dR}}(M)$ 或简写为 $Z^k(M)$。
  $$
  Z^k(M) = \{ \omega \in \Omega^k(M) \mid d\omega = 0 \}
  $$

* **恰当形式 (Exact Forms):** 一个 k-形式 $\omega \in \Omega^k(M)$ 称为 **恰当的 (exact)**，如果存在一个 (k-1)-形式 $\eta \in \Omega^{k-1}(M)$，使得 $\omega = d\eta$。所有恰当 k-形式构成的向量空间是 $\operatorname{im}(d: \Omega^{k-1} \to \Omega^k)$，记为 $B^k_{\text{dR}}(M)$ 或简写为 $B^k(M)$。
  $$
  B^k(M) = \{ \omega \in \Omega^k(M) \mid \exists \eta \in \Omega^{k-1}(M), \omega = d\eta \}
  $$

* 关键点：由于 $d^2 = 0$，**每一个恰当形式都是闭的**：$B^k(M) \subseteq Z^k(M)$。然而，**并非每一个闭形式都是恰当的**！闭形式是否恰当，取决于流形的拓扑结构。例如，在 $\mathbb{R}^2 \setminus \{(0,0)\}$ 上，形式 $\omega = \frac{-y dx + x dy}{x^2 + y^2}$ 是闭的 ($d\omega=0$)，但不是恰当的（如果它是恰当的，围绕原点的环路积分会是 0，但实际上是 $2\pi$）。

* **de Rham 上同调群 (de Rham Cohomology Groups):** 第 k 个 de Rham 上同调群 $H^k_{\text{dR}}(M)$ 定义为闭 k-形式模掉恰当 k-形式的商空间：
  $$
  H^k_{\text{dR}}(M) := Z^k(M) / B^k(M) = \frac{\ker(d: \Omega^k \to \Omega^{k+1})}{\operatorname{im}(d: \Omega^{k-1} \to \Omega^k)}
  $$

  *   元素是 **上同调类 (cohomology classes)** $[\omega]$，其中 $\omega$ 是一个闭形式（代表元）。两个闭形式 $\omega_1$ 和 $\omega_2$ 属于同一个上同调类（即 $[\omega_1] = [\omega_2]$）当且仅当它们的差是恰当的：$\omega_1 - \omega_2 = d\eta$（对于某个 $\eta \in \Omega^{k-1}(M)$）。
  *   $H^k_{\text{dR}}(M)$ 是一个实向量空间（因为形式系数在 $\mathbb{R}$ 上）。
  *   当 $k=0$ 时，$H^0_{\text{dR}}(M) = \ker(d: \Omega^0 \to \Omega^1)$。由于 $df=0$ 意味着 $f$ 在连通分支上是常数，所以 $H^0_{\text{dR}}(M) \cong \mathbb{R}^c$，其中 $c$ 是 $M$ 的连通分支数。
  *   当 $k > \dim M$ 时，$H^k_{\text{dR}}(M) = 0$。



### 2.1 看见”洞“：常见微分流形上的 de Rham 上同调群的计算

我们来计算 $ \mathbb{R}^n $、$ S^1 $、$ T^2 $ 的 de Rham 上同调 $ H^k_{dR} $，并通过微分形式与积分直观展示“洞”的拓扑信息。

#### 2.1.1. $ \mathbb{R}^n $ 的 de Rham 上同调

**Poincaré 引理**：对 $ \mathbb{R}^n $（或任何可缩空间），有：
$$
H^k_{dR}(\mathbb{R}^n) \cong 
\begin{cases} 
\mathbb{R}, & k = 0, \\
0, & k \geq 1.
\end{cases}
$$

- **$ k = 0 $**：$ \mathbb{R}^n $ 连通，故 $ H^0 \cong \mathbb{R} $（常函数）。
- **$ k \geq 1 $**：
  - 设 $ \omega $ 是闭形式（$ d\omega = 0 $），需证明 $ \omega $ 恰当（$ \omega = d\eta $）。
  - 构造同伦算子 $ H: \Omega^k \to \Omega^{k-1} $ 满足 $ \omega = d(H\omega) + H(d\omega) $（对闭形式即 $ \omega = d(H\omega) $）。
  - 例如在 $ \mathbb{R}^n $ 上，取径向同伦 $ H\omega = \int_0^1 t^{k-1} \iota_{\vec{r}} \omega(t\vec{r}) \, dt $，其中 $ \vec{r} = \sum x^i \partial_i $。
  - 因此所有闭形式均为恰当，$ H^k = 0 $（$ k \geq 1 $）。

**几何意义**：  
$ \mathbb{R}^n $ 无“洞”，故高阶上同调平凡。闭形式局部总可积（如 $ \omega = df $）。

#### **2.1.2. 圆周 $ S^1 $ 的 de Rham 上同调**

$$
H^k_{dR}(S^1) \cong 
\begin{cases} 
\mathbb{R}, & k = 0, 1, \\
0, & k \geq 2.
\end{cases}
$$

**计算**：

1. **$ k = 0 $**：同 $ \mathbb{R}^n $，$ H^0 \cong \mathbb{R} $（连通）。
2. **$ k = 1 $**：
   - 1-形式 $ \omega = g(\theta)d\theta $ 闭（因 $ d\omega = 0 $）。
   - $ \omega $ 恰当当且仅当 $ \int_{S^1} \omega = 0 $（即 $ g = df/d\theta $，则 $ \int g \, d\theta = f(2\pi) - f(0) = 0 $）。
   - 故 $ H^1 \cong \mathbb{R} $，生成元为 $ [d\theta] $，其积分 $ \int_{S^1} d\theta = 2\pi \neq 0 $。
3. **$ k \geq 2 $**：平凡（因 $ \dim S^1 = 1 $）。

**“洞”的体现**：  

- $ H^1 \neq 0 $ 对应 $ S^1 $ 的 **1 维孔洞**（非收缩环路）。  
- 生成元 $ d\theta $ 积分非零，反映绕圆周的全局非平凡性。

#### **2.1.3. 环面 $ T^2 = S^1 \times S^1 $ 的 de Rham 上同调**

$$
H^k_{dR}(T^2) \cong 
\begin{cases} 
\mathbb{R}, & k = 0, 2, \\
\mathbb{R}^2, & k = 1, \\
0, & k \geq 3.
\end{cases}
$$

**计算**（利用 Künneth 公式 $ H^\bullet(M \times N) \cong H^\bullet(M) \otimes H^\bullet(N) $）：

1. **$ k = 0 $**：$ T^2 $ 连通，$ H^0 \cong \mathbb{R} $。
2. **$ k = 1 $**：
   - $ H^1_{dR}(S^1) \cong \mathbb{R} $ 的生成元 $ d\theta_1 $、$ d\theta_2 $（两圆周坐标）。
   - $ H^1(T^2) \cong \mathbb{R}^2 $，基为 $ [d\theta_1], [d\theta_2] $。
   - 积分 $ \int_{\gamma_i} d\theta_j = 2\pi \delta_{ij} $（$ \gamma_i $ 为第 $ i $ 方向生成环路）。
3. **$ k = 2 $**：
   - 体积形式 $ d\theta_1 \wedge d\theta_2 $ 闭且非恰当（$ \int_{T^2} d\theta_1 \wedge d\theta_2 = (2\pi)^2 \neq 0 $）。
   - $ H^2 \cong \mathbb{R} $，生成元 $ [d\theta_1 \wedge d\theta_2] $。
4. **$ k \geq 3 $**：平凡（因 $ \dim T^2 = 2 $）。

**“洞”的体现**：

- **$ H^1 \cong \mathbb{R}^2 $**：对应两个独立非收缩环路（如经纬圈 $ \gamma_1, \gamma_2 $）。
- **$ H^2 \cong \mathbb{R} $**：反映 $ T^2 $ 的 **2 维空洞**（内部中空）。

#### **2.1.4. 总结对比**

| 流形             | $ H^0 $        | $ H^1 $          | $ H^2 $        | 拓扑解释                      |
| ---------------- | -------------- | ---------------- | -------------- | ----------------------------- |
| $ \mathbb{R}^n $ | $ \mathbb{R} $ | 0                | 0              | 无洞                          |
| $ S^1 $          | $ \mathbb{R} $ | $ \mathbb{R} $   | 0              | 1 个 1 维孔洞（环路）         |
| $ T^2 $          | $ \mathbb{R} $ | $ \mathbb{R}^2 $ | $ \mathbb{R} $ | 2 个 1 维孔洞 + 1 个 2 维空洞 |

**关键结论**：

- **$ H^1 $ 的维数** = 流形上“独立非收缩环路”数。  
- **$ H^2 $ 非零** = 流形有“封闭的内部空洞”（如 $ T^2 $ 的中空部分）。  
- **积分检测**：闭形式的积分是否为零判断其是否恰当（如 $ \int_{S^1} d\theta $、$ \int_{T^2} d\theta_1 \wedge d\theta_2 $）。  

这些计算清晰地展示了 de Rham 上同调如何通过微分形式捕捉流形的拓扑“孔洞”。



### 2.2 de Rham 上同调群的生成元

#### **2.2.1. 生成元的概念与几何意义**

在 de Rham 上同调群 $ H^k_{\text{dR}}(M) $ 中，**生成元**是代表上同调类的一组闭微分形式，它们满足：

- **闭形式**：$ d\omega = 0 $（局部可积条件）。
- **非恰当形式**：不存在 $ \eta $ 使得 $ \omega = d\eta $（全局非平凡性）。
- **线性无关**：它们的等价类在商空间中线性无关。
- **生成性**：任意闭形式的上同调类可表示为生成元的线性组合（模恰当形式）。

**几何意义**：  
生成元对应流形的**拓扑特征**（如“洞”或“空腔”）：

- $ H^1 $ 的生成元描述“1 维洞”（不可收缩的闭曲线）。
- $ H^2 $ 的生成元描述“2 维空腔”（不可收缩的闭曲面）。
  通过积分可量化这些特征（如 $ \int_{S^1} d\theta = 2\pi $）。

#### 2.2.**2. 为什么微分形式是生成元？**

微分形式作为生成元，本质是因其**积分值**和**拓扑不变性**反映了流形的全局结构。

##### **例 1：$ S^1 $ 的生成元 $ d\theta $**

- **闭形式**：$ d(d\theta) = 0 $（外微分的性质）。
- **非恰当**：$ \int_{S^1} d\theta = 2\pi \neq 0 $，但 Stokes 定理要求恰当形式在闭链上积分为零。
- **生成性**：任何闭 1-形式 $ \omega $ 满足 $ [\omega] = c [d\theta] $，其中 $ c = \frac{1}{2\pi} \int_{S^1} \omega $.

##### **例 2：$ S^2 $ 的生成元 $ \sin\phi \, d\theta \wedge d\phi $**

- **闭形式**：在 $ S^2 $ 上 $ d(\sin\phi \, d\theta \wedge d\phi) = 0 $（体积形式闭）。
- **非恰当**：$ \int_{S^2} \sin\phi \, d\theta \wedge d\phi = 4\pi \neq 0 $。
- **生成性**：由 Hodge 理论，$ H^2(S^2) $ 由调和形式生成，而球面面积形式是唯一的调和 2-形式。

##### **关键原因**：

1. **积分不变量**：生成元在闭链上的积分非零且拓扑不变（如 $ \int_{S^1} d\theta $ 不依赖坐标）。
2. **de Rham 定理**：$ H^k_{\text{dR}}(M) \cong H^k_{\text{sing}}(M; \mathbb{R}) $，生成元对应奇异上同调的生成元。
3. **几何实现**：生成元可由流形的对称性构造（如 $ S^1 $ 的旋转对称性给出 $ d\theta $，$ S^2 $ 的球坐标给出面积形式）。

#### 2.2.**3. 为什么 $ H^k_{\text{dR}}(M) $ 是一个群（而不仅是线性空间）？**

de Rham 上同调群同时具备**线性空间结构**和**群结构**：

- **向量空间结构**：它是商空间 $ \ker d / \operatorname{im} d $，其中：
  - $ \ker d $（闭形式）是实向量空间。
  - $ \operatorname{im} d $（恰当形式）是其子空间。
    因此 $ H^k_{\text{dR}}(M) $ 是实向量空间（线性空间）。

- **加法群结构**：作为向量空间，它自然是一个**加法阿贝尔群**：
  - **封闭性**：若 $ [\omega_1], [\omega_2] \in H^k_{\text{dR}}(M) $，则 $ [\omega_1] + [\omega_2] = [\omega_1 + \omega_2] $。
  - **结合律**：$ ([\omega_1] + [\omega_2]) + [\omega_3] = [\omega_1] + ([\omega_2] + [\omega_3]) $。
  - **单位元**：$ [0] $（恰当形式类）。
  - **逆元**：$ -[\omega] = [-\omega] $。

##### **为什么强调“群”而非仅“线性空间”？**

1. **历史与推广**：  
   - 上同调理论起源于代数拓扑中的群结构（如奇异上同调群 $ H^k_{\text{sing}}(M; G) $ 可定义在任意交换群 $ G $ 上）。
   - de Rham 上同调是特例（$ G = \mathbb{R} $），但结构定理表明它是**群与线性空间的统一体**。

2. **范畴观点**：  
   - 在范畴论中，上同调函子取值在**阿贝尔群范畴**中。
   - de Rham 上同调通过积分与奇异上同调的同构 $ H^k_{\text{dR}}(M) \cong \operatorname{Hom}(H_k(M), \mathbb{R}) $ 继承了群结构。

3. **物理与几何应用**：  
   - 群结构允许定义**上积**（cup product）$ \smile: H^k \times H^m \to H^{k+m} $，赋予 $ H^*_{\text{dR}}(M) $ 环结构。
   - 例如在 $ S^2 $ 上，$ H^0 \smile H^2 $ 生成 $ H^2 $，对应体积形式。

#### 2.2.4. **生成元与拓扑特征的对应表**

| 流形   | $ H^k $ | 生成元                               | 积分值                                             | 拓扑特征               |
| ------ | ------- | ------------------------------------ | -------------------------------------------------- | ---------------------- |
| $ S^1$ | $ H^0 $ | 常值函数 $ 1 $                       | $ \int_{S^1} 1 = 1 $                               | 连通性（0 维洞）       |
|        | $ H^1 $ | $ d\theta $                          | $ \int_{S^1} d\theta = 2\pi $                      | 1 维洞（不可收缩的圆） |
| $ T^2$ | $ H^1 $ | $ d\theta_1, d\theta_2 $             | $ \int_{C_i} d\theta_j = 2\pi \delta_{ij} $        | 两个独立的 1 维洞      |
|        | $ H^2 $ | $ d\theta_1 \wedge d\theta_2 $       | $ \int_{T^2} d\theta_1 \wedge d\theta_2 = 4\pi^2 $ | 2 维空腔（环面内部）   |
| $ S^2$ | $ H^2 $ | $ \sin\phi \, d\theta \wedge d\phi $ | $ \int_{S^2} = 4\pi $                              | 2 维空腔（球体内部）   |



### 2.3 Poincaré 引理与 Stokes 公式的关系及其揭示的“洞”结构

#### 2.3.**1. Poincaré 引理：局部无洞性**

设 $ U \subset \mathbb{R}^n $ 是**可收缩开集**（如星形区域），则：

- **内容**：对任意闭形式 $ \omega $（即 $ d\omega = 0 $)，存在形式 $ \eta $ 使得 $ \omega = d\eta $。

- **数学表述**：
  $$
  H^k_{\text{dR}}(U) = 0 \quad \text{当} \quad k \geq 1.
  $$

- **几何意义**：在局部可收缩区域中，闭形式必是恰当形式（即存在“原函数”），这反映了**局部无洞**的特性。

> **例**：在 $ \mathbb{R}^2 $ 中，闭形式 $ \omega = -y dx + x dy $ 满足 $ d\omega = 0 $。  
> Poincaré 引理保证存在函数 $ f = \frac{1}{2}(x^2 + y^2) $ 使得 $ df = -y dx + x dy $（星形区域上成立）。

#### 2.3.**2. Stokes 公式：全局积分与边界**

设 $ M $ 是 $ n $-维定向流形，$ \omega $ 是紧支集 $ (n-1) $-形式，则：
$$
\int_M d\omega = \int_{\partial M} \omega.
$$

- **几何意义**：流形内部的微分运算（外微分）在积分意义上等于边界上的通量。
- **核心**：建立了**内部微分**与**边界积分**的桥梁。

#### 2.3.**3. 二者关系：局部与全局的桥梁**

| **Poincaré 引理**      | **Stokes 公式**        | **关系**                      |
| ---------------------- | ---------------------- | ----------------------------- |
| 局部性质（可收缩开集） | 全局性质（带边流形）   | Stokes 是 Poincaré 的积分版本 |
| 断言闭形式局部可积     | 计算闭形式在闭链上积分 | 共同揭示拓扑障碍              |

##### **关键联系**：

1. **Poincaré 引理的证明依赖 Stokes 公式**：  
   构造同伦算子 $ K: \Omega^k(U) \to \Omega^{k-1}(U) $ 使得：
   $$
   \omega = d(K\omega) + K(d\omega)
   $$
   当 $ d\omega = 0 $ 时，$ \omega = d(K\omega) $。  
   **构造核心**：对线性路径积分，本质是应用 Stokes 公式于单形链。

2. **Stokes 公式解释 Poincaré 引理的失效**：  
   若闭形式 $ \omega $ 在闭链 $ c $ 上积分非零：
   $$
   \int_c \omega \neq 0, \quad \partial c = \emptyset
   $$
   则 $ \omega $ 不可能恰当（否则由 Stokes 公式，积分必为零）。  
   **失效原因**：流形存在“洞”，导致局部可积性被破坏。

#### 2.3.**4. 为什么这揭示了“洞”？**

##### **机制分析**：

1. **无洞情形（Poincaré 引理成立）**：  

   - 可收缩区域中，任意闭链 $ c $ 是某区域的边界（$ c = \partial D $）。  

   - 由 Stokes 公式：
     $$
     \int_c \omega = \int_{\partial D} \omega = \int_D d\omega = 0 \quad (\text{因 } d\omega=0).
     $$
     积分恒为零 → 闭形式必恰当。

2. **有洞情形（Poincaré 引理失效）**：  

   - 存在闭链 $ c $ 不能表示为边界（如 $ S^1 $ 中的圆周）。  

   - 若存在闭形式 $ \omega $ 满足：
     $$
     \int_c \omega \neq 0
     $$
     则 $ \omega $ 非恰当（否则 Stokes 公式强制积分为零）。  

   - **物理意义**：  
     “洞”阻碍了闭链收缩到点，导致积分路径不可忽略（如磁场沿环面的积分）。

##### **示例：圆周 $ S^1 $ 的洞**

- **闭形式**：$ \omega = d\theta $（满足 $ d\omega = 0 $)，但非全局恰当（$ \theta $ 非单值函数）。  

- **积分检测**：取闭链 $ c = S^1 $（整个圆周），
  $$
  \int_{S^1} d\theta = 2\pi \neq 0.
  $$

- **Stokes 公式分析**：  
  若 $ \omega $ 恰当（即 $ \omega = df $)，则：
  $$
  \int_{S^1} df = f(2\pi) - f(0) = 0 \quad (\text{矛盾}).
  $$
  故 $ \omega $ 非恰当 → 反映 $ S^1 $ 有 1 维洞。

##### **示例：球面 $ S^2 $ 的空腔**

- **闭形式**：体积形式 $ \omega = \sin\phi  d\phi \wedge d\theta $（满足 $ d\omega = 0 $)。  

- **积分检测**：取闭链 $ c = S^2 $（整个球面），
  $$
  \int_{S^2} \omega = 4\pi \neq 0.
  $$

- **结论**：$ \omega $ 非恰当 → 反映 $ S^2 $ 有 2 维空腔。

#### 2.3.5 总结：Poincaré-Stokes 框架下的“洞”理论

| **概念**          | **数学描述**                                  | **揭示的拓扑结构**      |
| ----------------- | --------------------------------------------- | ----------------------- |
| **Poincaré 引理** | $ H^k_{\text{dR}}(U) = 0 $                    | 局部无洞（可收缩区域）  |
| **Stokes 公式**   | $ \int_M d\omega = \int_{\partial M} \omega $ | 全局积分与边界的关系    |
| **洞的检测**      | $ \exists  c: \int_c \omega \neq 0 $          | 闭链非边界 → $ k $-维洞 |

**核心结论**：  

- **Poincaré 引理** 表明：局部无洞时，闭形式必可积（恰当）。  
- **Stokes 公式** 提供检测工具：当闭形式沿闭链积分非零时，说明该闭链无法收缩到点（即存在“洞”），导致 Poincaré 引理失效。  
- **de Rham 上同调** $ \dim H^k_{\text{dR}}(M) $ 直接给出 $ k $-维洞的数量（Betti 数）。  

这一框架将局部微分性质（Poincaré 引理）与全局积分不变量（Stokes 公式）统一，揭示了微分形式如何捕捉流形的拓扑障碍（“洞”）。



### 2.4. de Rham 上同调的几何与拓扑意义

* **测量“洞” (Measuring "Holes"):** de Rham 上同调群的维度（称为 **Betti 数 (Betti numbers)** $b_k := \dim H^k_{\text{dR}}(M)$）提供了流形上不同维数的“洞”的信息。

  *   $b_0$：连通分支数（0维洞）。
  *   $b_1$：本质上的“一维洞”的数量（例如，圆圈有一个，环面有两个）。
  *   $b_2$：本质上的“封闭曲面”或“空腔”的数量（例如，球面有一个，环面有一个）。
  *   以此类推。

* **上同调类的不变性 (Invariance of Cohomology Classes):** 关键定理表明，de Rham 上同调群是 **微分同胚不变量 (diffeomorphism invariant)**，甚至是 **同伦不变量 (homotopy invariant)**。这意味着：

  1.  如果两个流形 $M$ 和 $N$ 是 **微分同胚 (diffeomorphic)** 的（即存在光滑的双射，其逆也是光滑的），那么 $H^k_{\text{dR}}(M) \cong H^k_{\text{dR}}(N)$ 对所有 $k$ 成立。
  2.  更一般地，如果两个流形是 **同伦等价 (homotopy equivalent)** 的（即一个可以连续形变到另一个），那么它们的 de Rham 上同调群也是同构的。

  *   因此，de Rham 上同调群是流形的 **光滑同伦型 (smooth homotopy type)** 的不变量，是强大的拓扑不变量。

* **de Rham 定理 (de Rham's Theorem):** 这是整个理论的巅峰之作。它建立了光滑微分形式构造的 de Rham 上同调与通过单纯复形或 Čech 上同调等纯拓扑方法定义的 **奇异上同调 (singular cohomology)** 之间的深刻联系：
  $$
  \boxed{H^k_{\text{dR}}(M; \mathbb{R}) \cong H^k_{\text{sing}}(M; \mathbb{R})}
  $$

  *   这个定理断言，de Rham 上同调群（用微分形式定义）与系数在实数 $\mathbb{R}$ 上的奇异上同调群是同构的。
  *   意义：它证明了 de Rham 上同调确实捕获了流形的 **真实拓扑本质 (intrinsic topological essence)**，而不仅仅是光滑结构的信息。光滑微分形式提供了一种强大的 **计算方法 (computational tool)** 来访问和计算拓扑不变量（奇异上同调群）。

* **Stokes 定理的推广 (Generalization of Stokes' Theorem):** 经典的 Stokes 定理（格林定理、高斯散度定理、标准 Stokes 定理的统一）可以优雅地表述为：
  $$
  \int_M d\omega = \int_{\partial M} \omega
  $$
  其中 $M$ 是带边界的紧定向流形，$\dim \omega = \dim M - 1$，$\partial M$ 是其边界。这个定理深刻地将外微分 $d$（局部算子）与积分（全局算子）联系起来，并且是证明 de Rham 定理的关键工具之一。它体现了局部微积分操作与全局拓扑性质之间的深刻联系。

* **函子性 (Functoriality):** 光滑映射 $f: M \to N$ 通过 **拉回 (pullback)** $f^*: \Omega^k(N) \to \Omega^k(M)$ 诱导线性映射 $f^*: H^k_{\text{dR}}(N) \to H^k_{\text{dR}}(M)$。拉回与外微分交换 ($f^* \circ d_N = d_M \circ f^*$)，所以它将闭形式映到闭形式，将恰当形式映到恰当形式，从而诱导上同调映射。这使得 $H^*_{\text{dR}}$ 成为一个 **反变函子 (contravariant functor)**。

* **同伦不变性 (Homotopy Invariance):** 如果两个光滑映射 $f, g: M \to N$ 是 **光滑同伦 (smoothly homotopic)** 的（即存在光滑映射 $H: M \times [0,1] \to N$ 使得 $H(x,0)=f(x)$, $H(x,1)=g(x)$），则它们诱导相同的上同调映射：$f^* = g^*: H^k_{\text{dR}}(N) \to H^k_{\text{dR}}(M)$。这直接推出同伦等价流形有同构的上同调。

de Rham 上同调是微分几何和拓扑学中一座辉煌的桥梁：

1.  **局部到全局 (Local to Global):** 它利用定义在流形局部坐标卡上的光滑微分形式（及其外微分算子 $d$）构造出刻画流形整体拓扑结构的向量空间 $H^k_{\text{dR}}(M)$。
2.  **分析到拓扑 (Analysis to Topology):** 它通过分析工具（微分形式、积分、微分方程）计算和揭示了纯粹的拓扑不变量（Betti 数、上同调群）。
3.  **de Rham 定理的核心地位:** 该定理确立了这种构造的有效性，证明它与任何其他方法定义的实系数上同调理论（如奇异上同调）本质上是相同的。
4.  **强大的工具 (Powerful Tool):** 它提供了研究流形拓扑的强有力工具（函子性、同伦不变性、Mayer-Vietoris 序列、Hodge 理论等），在数学物理（如规范理论、量子场论）、复几何、辛几何等领域有广泛应用。
5.  **几何内涵 (Geometric Insight):** 上同调类代表了存在“拓扑障碍”阻止闭形式成为恰当形式。这种障碍对应于流形中不同维数的“洞”或“不可收缩的圈/曲面”。积分（Stokes 定理）是探测这些类的重要方式。

总而言之，de Rham 上同调深刻揭示了光滑流形上局部微分结构与全局拓扑性质之间密不可分的联系，是理解现代几何与拓扑不可或缺的语言和工具。



## 3. 同伦群与同调群

### 3.1. 同伦群

**一、同伦的基本概念**

1. **同伦（Homotopy）的定义**
   设 $X$ 为拓扑空间，$\;f,g: A\to X$ 为两个连续映射，其中 $A$ 也是拓扑空间。若存在连续映射

   $$
     H: A\times[0,1]\;\longrightarrow\;X
   $$

   使得

   $$
     H(a,0)=f(a),\quad H(a,1)=g(a)\quad(\forall\,a\in A),
   $$

   则称 $H$ 是从 $f$ 到 $g$ 的**同伦**，记作 $f\simeq g$。

   * 直观理解：$H$ 把 $f$ 的图像在参数 $t$ （“时间”） 上连续变形到 $g$ 的图像。
   * 特殊情况：若 $A=[0,1]$ 且 $f,g$ 都是基点回路（即 $f(0)=f(1)=g(0)=g(1)=x_0$），则我们谈的是**基点回路**之间的同伦。

2. **同伦的等价关系性质**
   对于固定域 $A\to X$ 的所有连续映射，同伦关系 $\simeq$ 满足：

   * **自反性**：显然，可取 $H(a,t)=f(a)$（与时间无关），则 $f\simeq f$。

   * **对称性**：若 $H$ 是从 $f$ 到 $g$ 的同伦，则定义 $H'(a,t)=H(a,1-t)$，即可得到从 $g$ 到 $f$ 的同伦，故 $f\simeq g\implies g\simeq f$。

   * **传递性**：若 $H_1$ 是 $f\simeq g$ 的同伦，$H_2$ 是 $g\simeq h$ 的同伦，可将它们“拼接”构造出 $f\simeq h$：

     $$
       H(a,t)=
       \begin{cases}
         H_1(a,2t), & 0\le t\le \tfrac12,\\
         H_2(a,2t-1), & \tfrac12\le t\le1.
       \end{cases}
     $$

     这样 $H(a,0)=f(a)$，$H(a,1)=h(a)$，完成传递性。

   由此，同伦在映射集合上形成一个等价关系，可将映射按同伦类（equivalence class）来分类。

**二、第一同伦群 $\pi_1(X,x_0)$**

1. **基点回路与同伦类**
   选定基点 $x_0\in X$。考虑所有满足

   $$
     \gamma:\,[0,1]\;\to\;X,\quad \gamma(0)=\gamma(1)=x_0
   $$

   的闭路，将两条闭路若存在保持基点不动的同伦即视为同一类。同伦类的集合记为 $\pi_1(X,x_0)$。

2. **群运算：回路的拼接**
   对两条回路的同伦类 $[\alpha]$ 与 $[\beta]$，定义拼接 $\alpha*\beta$：

   $$
     (\alpha*\beta)(t)=
     \begin{cases}
       \alpha(2t), & 0\le t\le\tfrac12,\\
       \beta(2t-1),& \tfrac12\le t\le1.
     \end{cases}
   $$

   这一操作在同伦类意义下良定义，并满足群的四条公理：

   * **封闭性**：拼接仍是基点回路；
   * **结合律**：$(\alpha*\beta)*\gamma\simeq \alpha*(\beta*\gamma)$（可构造自然的三段时间拼接同伦）；
   * **单位元**：恒等回路 $e(t)=x_0$ 的同伦类；
   * **逆元**：回路 $\alpha$ 的反向 $\alpha^{-1}(t)=\alpha(1-t)$。

   因此，$\pi_1(X,x_0)$ 构成一个群，称为**第一同伦群**或**基本群**（Fundamental Group）。

**三、高阶同伦群 $\pi_n(X,x_0)$**

1. **定义**
   对 $n\ge2$，将区间 $[0,1]$ 替换为 $n$-维立方体 $[0,1]^n$，并把整个边界 $\partial[0,1]^n$ 均映到基点 $x_0$。具体地，考察连续映射

   $$
     f:\;([0,1]^n,\,\partial[0,1]^n)\;\longrightarrow\;(X,\;x_0),
   $$

   即

   $$
     f(\,x\,)=x_0,\quad\forall\,x\in\partial[0,1]^n.
   $$

   若两个此类映射可通过保持边界不动的同伦连通，则属于同一同伦类，记为 $\pi_n(X,x_0)$。

2. **群结构**
   当 $n\ge2$ 时，$\pi_n(X,x_0)$ 在同伦类层面也可以用“拼接”定义加法：

   * 将立方体的一个半（比如第一坐标 $x_1\in[0,\tfrac12]$）用于 $f$，另一半用于 $g$，构造拼接映射。
   * 对同伦类，这一操作满足阿贝尔群（交换性）的所有性质。

3. **代替定义：基于球面**
   也可等价地将 $[0,1]^n$ 及其边界压扁为 $n$-维球面 $S^n$，考察基点映射

   $$
     f:(S^n,*)\to (X,x_0)
   $$

   的同伦类，得到相同的 $\pi_n(X,x_0)$。

**四、典型例子**

1. **欧几里得空间 $\mathbb{R}^m$**

   $$
     \pi_n(\mathbb{R}^m)=0,\quad\forall\,n\ge1.
   $$

   因为 $\mathbb{R}^m$ 为可收缩空间，可同伦到一个点，所有回路和球面映射都可收缩。

2. **圆环 $S^1$**

   * $\pi_1(S^1)\cong\mathbb{Z}$：每个整数对应绕圆 “转”的圈数与方向；
   * $\pi_n(S^1)=0$ 对所有 $n\ge2$。

3. **高维球面 $S^n$**

   * $\pi_n(S^n)\cong\mathbb{Z}$：把 $S^n$ 自身映射到 $S^n$ 的“度数”（degree）；
   * $\pi_k(S^n)=0$ 当 $k<n$。
   * 但当 $k>n$ 时，同伦群变得极其复杂（例如 $\pi_3(S^2)\cong\mathbb{Z}$，对应 Hopf 纤维映射）。

4. **环面 $T^2=S^1\times S^1$**

   * $\pi_1(T^2)\cong\pi_1(S^1)\times\pi_1(S^1)\cong\mathbb{Z}\times\mathbb{Z}$；
   * $\pi_n(T^2)=0\,(n\ge2)$。



### 3.2 同调群

**一、单形 (Simplex) 与链 (Chains)**

1. **单形的定义**

   * 一个 **$n$-单形** $\sigma$ 是 $n+1$ 个**仿射无关**点（顶点）$\{v_0,v_1,\dots,v_n\}$ 在欧氏空间中的凸包：

     $$
       \sigma=[v_0,v_1,\dots,v_n]
       =\Bigl\{\sum_{i=0}^n t_i\,v_i\;\Big|\;t_i\ge0,\;\sum_{i=0}^n t_i=1\Bigr\}.
     $$

   * 特殊情形：

     * $0$-单形就是一个点 $[v_0]$。
     * $1$-单形是线段 $[v_0,v_1]$。
     * $2$-单形是三角形 $[v_0,v_1,v_2]$。

2. **链群 $C_n(X)$**

   * 给定一个有向三角剖分（或更一般的**单形复形**）$X$，令所有 $n$-单形的形式线性组合构成自由阿贝尔群

     $$
       C_n(X)=\Bigl\{\sum_k a_k\,\sigma_k\;\Big|\;a_k\in\mathbb{Z},\;\sigma_k\text{ 是 }n\text{-单形}\Bigr\}.
     $$

   * 元素称为**$n$链**。

   * 这里系数群选 $\mathbb{Z}$；也可选用其他系数环（如 $\mathbb{Z}_2$、$\mathbb{R}$ 等）。

**二、边界算子 (Boundary Operator)**

1. **定义**
   对于一个 $n$-单形 $\sigma=[v_0,v_1,\dots,v_n]$，其 **边界** $\partial_n\sigma$ 定义为：

   $$
     \partial_n[v_0,v_1,\dots,v_n]
     =\sum_{i=0}^n(-1)^i\,[v_0,\dots,\widehat{v_i},\dots,v_n],
   $$

   其中 $\widehat{v_i}$ 表示去掉顶点 $v_i$，得到一个 $(n-1)$-单形。

2. **对链的延拓**
   $\partial_n$ 在线性意义下作用于整个链群：

   $$
     \partial_n\Bigl(\sum_k a_k\,\sigma_k\Bigr)
     =\sum_k a_k\,\partial_n\sigma_k.
   $$

3. **性质：$\partial_{n-1}\circ\partial_n=0$**

   * 计算可验证对任何单形
     $\partial_{n-1}(\partial_n\sigma)=0$。
   * 直观上，“边界的边界为空”：三角形的边界是三条有向线段，这三条线段拼起来是一个闭环，同其边界（端点）互相抵消。

**三、链复形 (Chain Complex)**

1. **定义**
   将各维链群及其边界算子连接起来，构成一个复形：

   $$
     \cdots 
     \;\xrightarrow{\;\partial_{n+1}\;}C_n(X)
     \;\xrightarrow{\;\partial_n\;}C_{n-1}(X)
     \;\xrightarrow{\;\partial_{n-1}\;} 
     \cdots 
     \;\xrightarrow{\;\partial_1\;}C_0(X)
     \;\xrightarrow{\;\partial_0\;}0.
   $$

   满足 $\partial_{n-1}\circ\partial_n=0$。

2. **子群：循环与边界**

   * **$n$-循环群**

     $$
       Z_n(X)=\ker(\partial_n)
       =\{\,c\in C_n(X)\mid \partial_n c=0\}.
     $$

     元素称为“无边界”的 $n$-链。

   * **$n$-边界群**

     $$
       B_n(X)=\mathrm{im}(\partial_{n+1})
       =\{\,\partial_{n+1}c\mid c\in C_{n+1}(X)\}.
     $$

     元素称为“被更高维单形填充”的 $n$-链。

   * 由于 $\partial_{n}\circ\partial_{n+1}=0$，总有 $B_n(X)\subseteq Z_n(X)$。

**四、同调群 (Homology Groups) 的定义**

对每一个维度 $n$，定义**第 $n$ 同调群**为商群：

$$
  H_n(X)
  =\frac{Z_n(X)}{B_n(X)}.
$$

它度量了“无法用更高维单形填充的”$n$-洞：

* **分子** $Z_n$：所有“环”（无边界）
* **分母** $B_n$：所有真的是更高维边界的环

**五、典型例子**

1. **圆环 $S^1$**

   * $C_0$ 由顶点生成，$C_1$ 由有向边段生成。

   * 计算可得

     $$
       H_0(S^1)\cong\mathbb{Z},\quad
       H_1(S^1)\cong\mathbb{Z},\quad
       H_n(S^1)=0\;(n\ge2).
     $$

   * 其中 $H_1\cong\mathbb{Z}$ 对应“绕圈一次”的类；$H_0\cong\mathbb{Z}$ 对应连通分支数。

2. **二维球面 $S^2$**

   * 三角剖分后可求得

     $$
       H_0(S^2)\cong\mathbb{Z},\quad
       H_2(S^2)\cong\mathbb{Z},\quad
       H_1(S^2)=0,\quad H_n(S^2)=0\;(n\ge3).
     $$

   * $H_2$ 对应球面的“面洞”；没有一维洞。

3. **环面 $T^2=S^1\times S^1$**

   * 可用积空间的 Künneth 公式（或直接剖分）算得

     $$
       H_0(T^2)\cong\mathbb{Z},\quad
       H_1(T^2)\cong\mathbb{Z}^2,\quad
       H_2(T^2)\cong\mathbb{Z},\quad
       H_n(T^2)=0\;(n\ge3).
     $$

   * $H_1\cong\mathbb{Z}^2$ 分别对应环面上“经圈”“纬圈”两种不可填充的回路；$H_2$ 对应整体的面洞。

**六、小结**

* **单形** 是构造空间的基本“砖块”；
* **边界算子** 给出砖块之间的拼接关系，并满足“边界的边界为空”；
* **链复形** 把所有维度串连起来，通过核与像的商而得同调群；
* **同调群** $H_n(X)$ 则精确刻画了空间中**$n$ 维洞**的“数量”与“类型”。

掌握这一套从几何分解到代数商群的机制，是代数拓扑的核心工具，也是后续研究胞腔复形、相对同调、持续同调等更高阶理论的基础。



### 3.3 同伦群、同调群、de Rham 上同调群之间的关系

下面从几个层面来把三者串联起来，重点突出它们在拓扑不变量提取上的相互作用。

#### 3.3.1 从基本群到第一同调群：阿贝尔化

* **基本群** $\pi_1(X,x_0)$ 是非交换群，刻画了基点回路的拼接结构。

* **第一同调群** $H_1(X)$ 始终是阿贝尔群，其本质就是把 $\pi_1$ “强制交换”以后的结果：

  $$
    H_1(X)\;\cong\;\pi_1(X,x_0)\big/\big[\pi_1,\pi_1\big],
  $$

  其中 $[\pi_1,\pi_1]$ 是由所有交换子生成的子群。

* **直观**：在同调中，只关心回路“绕过洞”的次数，不再区分先绕哪条、后绕哪条，于是群运算变为加法，也即阿贝尔化。

#### 3.3.2. de Rham 上同调：微分形式的计算流形不变量

* 对光滑流形 $M$ 而言，**de Rham 复形**

  $$
    0\;\longrightarrow\;\Omega^0(M)
    \xrightarrow{d}\;\Omega^1(M)
    \xrightarrow{d}\;\cdots
    \xrightarrow{d}\;\Omega^n(M)
    \;\longrightarrow\;0
  $$

  及其商群
  $\displaystyle H^k_{\rm dR}(M)=\ker d\big/\mathrm{im}\,d$
  给出了微分形式层面的“上同调”。

* **de Rham 定理**：

  $$
  H^k_{\rm dR}(M)\;\cong\;H^k_{\rm sing}(M;\mathbb{R}),
  $$

  即解析（微分形式）与拓扑（奇异上同调）的完全等价。



### 3.4 欧拉示性数

**一、欧拉示性数的定义**
对于一个有限维链复形（或具有良好三角剖分的紧致空间）$X$，**欧拉示性数**$\chi(X)$ 定义为链群维数的交错和：

$$
  \chi(X)\;=\;\sum_{k\ge0}(-1)^k\,\dim C_k\;=\;\sum_{k\ge0}(-1)^k\,(\#\text{\(k\)-单形})\,.
$$

由于“边界的边界为零”性质，可证明这一数值与具体剖分无关，是拓扑不变量。

**二、在同调群中的体现**
利用链复形到同调群的商结构，欧拉示性数也可表示为同调群维数（即**Betti 数**）的交错和：
$$
  \chi(X)
  =\sum_{k\ge0}(-1)^k\,\dim H_k(X)\,,
$$

其中

$$
  H_k(X)=\frac{\ker(\partial_k)}{\mathrm{im}(\partial_{k+1})},
  \qquad b_k=\dim H_k(X)
$$

称为第 $k$ 个 Betti 数。

> **证明要点**：
>
> $$
> \dim C_k \;=\;\dim Z_k + \dim B_{k-1},
> \quad
> \dim Z_k \;=\;\dim H_k + \dim B_k,
> $$
>
> 代入交错和后，$\dim B_k$ 相邻项抵消，余下 $\dim H_k$ 的交错和。

* **示例**

  * 对于二维球面 $S^2$：

    $$
      (b_0,b_1,b_2)=(1,0,1)\quad\Longrightarrow\quad
      \chi(S^2)=1-0+1=2.
    $$

  * 对于环面 $T^2$：

    $$
      (b_0,b_1,b_2)=(1,2,1)\quad\Longrightarrow\quad
      \chi(T^2)=1-2+1=0.
    $$

**三、在 de Rham 上同调群中的体现**
对于光滑流形 $M$，**de Rham 上同调群**
$$
  H^k_{\rm dR}(M)
  =\frac{\ker\bigl(d:\Omega^k\to\Omega^{k+1}\bigr)}
         {\mathrm{im}\bigl(d:\Omega^{k-1}\to\Omega^k\bigr)}
$$

也是有限维实向量空间，其维数即 **de Rham Betti 数**

$$
  b^k_{\rm dR}=\dim H^k_{\rm dR}(M).
$$

因 **de Rham 定理** 有自然同构

$$
  H^k_{\rm dR}(M)\;\cong\;H^k_{\rm sing}(M;\mathbb{R})
  \quad\Longrightarrow\quad
  b^k_{\rm dR}=b_k,
$$

于是同样可定义

$$
  \chi(M)
  =\sum_{k\ge0}(-1)^k\,b^k_{\rm dR}
  =\sum_{k\ge0}(-1)^k\,b_k.
$$

* **直观说明**：de Rham 共homology 提供了微分形式的“洞”计数，与奇异同调上同调等价，因而其 Betti 数与同调群 Betti 数相同，交错和不变。

**四、小结与延伸**

1. **链群视角**：$\chi$ 是单形数的交错和；
2. **同调群视角**：$\chi$ 是同调群 Betti 数的交错和；
3. **de Rham 视角**：同理为 de Rham 上同调群 Betti 数的交错和；
4. **一致性来源**：链复形上的核像维数交错抵消；de Rham 定理保证同调与上同调 Betti 数一致。

通过这三种等价定义，欧拉示性数成为连接几何分解、代数同调与解析微分形式的核心桥梁。



## 4. 延伸

### 4.1. Gauss-Bonnet-Chern

下面先从二维曲面说起，再概述高维的 Gauss–Bonnet–Chern 广义。

#### 4.1.1. 二维曲面的 Gauss–Bonnet 公式

设 $M$ 是一条光滑、有向、紧致的二维曲面（带或不带边界），记其高斯曲率为 $K$，第一基本形式度量所给的面积元素为 $\mathrm{d}A$，边界 $\partial M$ 的向外单位法线确定的边界曲线的测地曲率为 $k_g$，弧长元为 $\mathrm{d}s$。

**无边界情况**：

若 $\partial M=\emptyset$，则

$$
  \boxed{\int_M K\,\mathrm{d}A \;=\; 2\pi\,\chi(M),}
$$

其中 $\chi(M)$ 是欧拉示性数。

* **几何意义**：曲面上的曲率“总和”正比于其拓扑类型；例如球面 $\chi=2$，环面 $\chi=0$。

**带边界情况**：

若 $\partial M\neq\emptyset$，则

$$
  \boxed{\int_M K\,\mathrm{d}A \;+\;\int_{\partial M}k_g\,\mathrm{d}s
    \;=\;2\pi\,\chi(M).}
$$

* $\displaystyle\int_{\partial M}k_g\,\mathrm{d}s$ 称为**边界项**，度量边界相对于自身测地线的“偏离”量。
* 当边界为一系列光滑封闭曲线时，上式同样成立。

#### 4.1.2. 从三角形缺角看公式

* **三角形缺角公式**：在高斯曲率恒定为 $K$ 的区域内，任意 geodesic 三角形的内角和满足

  $$
    \alpha+\beta+\gamma \;=\;\pi\;+\;\int_{\text{三角形}}K\,\mathrm{d}A.
  $$

* 将整个曲面划分为一系列 geodesic 多边形，累加缺角，就能得到全局的 Gauss–Bonnet。

#### 4.1.3. 微分形式与 Chern–Weil 理论视角

把二维曲面当作有向 Riemann 流形来看，其 Levi-Civita 连接给出局部正交框 $\{e_1,e_2\}$ 下的连接 1-形式 $\omega$，其曲率 2-形式为

$$
  \Omega \;=\; \mathrm{d}\omega.
$$

经典计算可得

$$
  \Omega \;=\; K\,\mathrm{d}A.
$$

Gauss–Bonnet 可写为

$$
  \int_M\Omega \;=\;2\pi\,\chi(M).
$$

这正是 Euler 整数类（Euler class）在同调群 $\;H^2(M)$ 中的积分表述。

#### 4.1.4. 高维 Gauss–Bonnet–Chern 定理

对一个 $2n$ 维紧致无边界黎曼流形 $(M^{2n},g)$，令其曲率张量对应的 2-形式矩阵为 $\{\Omega^i_j\}$，则其**Pfaffian** $\operatorname{Pf}(\tfrac{\Omega}{2\pi})$ 是顶维形式。Gauss–Bonnet–Chern 定理断言：

$$
\boxed{
    \int_{M^{2n}}\!\operatorname{Pf}\!\Bigl(\tfrac{\Omega}{2\pi}\Bigr)
    \;=\;\chi\bigl(M^{2n}\bigr).
  }
$$

* 当 $n=1$ 时，$\operatorname{Pf}(\tfrac{\Omega}{2\pi})=\tfrac{K}{2\pi}\mathrm{d}A$，恢复二维公式。
* 这是 Chern–Weil 理论中 Euler 类的具体化：Euler 类是一个位于 $H^{2n}(M)$ 的特征类，其上同调代表恰好给出 Pfaffian。

#### 4.1.5. 小结

1. **二维**：$\displaystyle\int_MK\,\mathrm{d}A\,+\!\int_{\partial M}k_g\,\mathrm{d}s=2\pi\chi(M)$。
2. **微分形式**：$\int_M\Omega=2\pi\chi(M)$，其中 $\Omega$ 是曲率 2-形式。
3. **高维**：$\displaystyle\int_M\operatorname{Pf}(\tfrac{\Omega}{2\pi})=\chi(M)$，关联 Euler 类。

Gauss–Bonnet（及其高维推广）将**几何**（曲率）与**拓扑**（Euler 特征）紧密联结，是全局微分几何的基石。



### 4.2 pushforward

在黎曼流形的语境中，“pushforward” 最常指一个光滑映射在切丛层面诱导出的微分（有时也称为 **导子** 或 **微分映射**），以及它在更高级结构（如度量、测度）上的推广。下面我们分层次、带例子地深入阐述。

#### 4.2.1. 定义：映射的微分（dφ 或 φ\_\*）

设

$$
  \phi\colon (M,g)\;\longrightarrow\;(N,h)
$$

是一光滑映射，其中 $(M,g)$、$(N,h)$ 是黎曼流形。对于任意点 $p\in M$，

* **切空间** $T_pM$ 是在 $p$ 处的所有速度向量的集合，

* **pushforward**（微分）记作

  $$
    \phi_{*\,p}\;=\;d\phi_p\;:\;T_pM\;\longrightarrow\;T_{\phi(p)}N,
  $$

  它将一个在 $M$ 上的切向量 $v\in T_pM$ “推”到 $N$ 上的切向量
  $\phi_{*\,p}(v)\in T_{\phi(p)}N$。

#### 4.2.2. 坐标表达

若在 $M$ 的邻域取局部坐标 $(x^1,\dots,x^m)$，在 $N$ 的邻域取坐标 $(y^1,\dots,y^n)$，且

$$
  \phi:\;x\longmapsto y=\bigl(y^i(x)\bigr)_{i=1}^n,
$$

则在点 $p$ 处的微分矩阵（Jacobian）为

$$
  \bigl(d\phi_p\bigr)^i{}_j
  =\frac{\partial y^i}{\partial x^j}\Big|_{p},
$$

对应于

$$
  \phi_{*\,p}\bigl(\partial_{x^j}\bigr)
  =\sum_{i=1}^n \frac{\partial y^i}{\partial x^j}\Big|_{p}\,\partial_{y^i}.
$$

#### 4.2.3. 基本性质

1. **线性映射**：对每个 $p$，$\phi_{*\,p}$ 是一个线性映射。

2. **复合映射所诱导的微分**：

   $$
     (\psi\circ\phi)_{*\,p}
     =\psi_{*\,\phi(p)}\;\circ\;\phi_{*\,p}.
   $$

3. **单位映射**：$\mathrm{id}_M{}_*{}_p$ 是恒等映射。

#### 4.2.4. 向量场与 pushforward

* 若 $X$ 是 $M$ 上的向量场（即 $p\mapsto X_p\in T_pM$），则

  $$
    (\phi_* X)_{\phi(p)} \;=\;\phi_{*\,p}(X_p)
  $$

  定义了一个 $N$ 上的向量场**前推**，但要保证其平滑性，通常需要 $\phi$ 是**微分同构**（局部双射）或至少是**亚浸（submersion）**。

* 若 $\phi$ 不是满射／单射，$\phi_*X$ 可能无法在全域定义为光滑且单值。

#### 4.2.5. 与 pullback 的对偶关系

在切空间与余切空间（cotangent space）之间，$\phi_*$ 与“pullback” $\phi^*$（作用于微分形式）彼此对偶：

$$
  \bigl(\phi^*\omega\bigr)_p(v)
  =\omega_{\phi(p)}\bigl(\phi_{*\,p}(v)\bigr),
  \quad
  \forall\,\omega\in T^*_{\phi(p)}N,\;v\in T_pM.
$$

这保证了微分形式与向量的配对在映射前后保持一致。

#### 4.2.6. 在黎曼度量下的特殊情形：等距与度量推送

若 $\phi$ 是**等距映射**（isometry），即

$$
  h_{\phi(p)}\bigl(\phi_{*\,p}(u),\,\phi_{*\,p}(v)\bigr)
  =g_p(u,v),
  \quad \forall\,u,v\in T_pM,
$$

则 $\phi_*$ 不仅保持切向量的光滑结构，还严格保持了“长度”和“夹角”。这是几何上最重要的 pushforward 情形。

#### 4.2.7. 推送 Riemannian 体积与测度

黎曼度量 $g$ 在 $M$ 上定义了体积形式 $\mathrm{d}V_g$。映射 $\phi$ 诱导出的 **pushforward 测度** $\phi_*(\mathrm{d}V_g)$ 定义为

$$
  \bigl(\phi_*\mathrm{d}V_g\bigr)(U)
  =\int_{\phi^{-1}(U)}\mathrm{d}V_g,
  \quad U\subset N
$$

是 $N$ 上的测度。在局部坐标下，其密度变化由 Jacobi 行列式
$\bigl|\det d\phi_p\bigr|$ 调整：
$$
  \phi_*\mathrm{d}V_g
  =\bigl|\det(d\phi)\bigr|\;\mathrm{d}V_h.
$$

#### 4.2.8. 小结

1. **pushforward** $\phi_{*\,p}=d\phi_p$ 将 $M$ 上的切向量线性推到 $N$ 上；
2. 它与 pullback $\phi^*$ 在代数上对偶；
3. 在黎曼几何中，等距映射的 pushforward 保持度量，测度的 pushforward 由 Jacobian 行列式刻画；
4. 各种几何工具——如第二基本形式、指数映射、Lie 群对称——都依赖对 pushforward 的深入理解。

掌握这些概念后，你就可以在流形映射、测地比较、子流形几何等多种场景中游刃有余地应用 pushforward。



### 4.3. pullback

在黎曼流形 $(M,g)$ 与 $(N,h)$ 之间，若有光滑映射

$$
  \phi: M \;\longrightarrow\; N,
$$

我们可以在多种几何对象上定义 **pullback**（拉回）操作，将 $N$ 上的对象“拉”回到 $M$ 上。下面分层次细述其定义、性质与典型应用。

#### 4.3.1. Pullback 的基本定义

1. **函数的 pullback**
   如 $f\colon N\to\Bbb R$ 是一个光滑函数，则

   $$
     \phi^*f \;=\; f\circ\phi
     \;:\;M\;\longrightarrow\;\Bbb R
   $$

   也是光滑函数。

2. **微分形式的 pullback**

   * **1-形式**：若 $\omega\in\Omega^1(N)$，在点 $p\in M$ 上定义
     $$
       (\phi^*\omega)_p(v)
       = \omega_{\phi(p)}\bigl(\phi_{*\,p}(v)\bigr),
       \qquad v\in T_pM.
     $$

     这样 $\phi^*\omega$ 是 $M$ 上的 1-形式。

   * **高阶形式**：对任意 $k$-形式 $\alpha\in\Omega^k(N)$，在 $p$ 处

     $$
       (\phi^*\alpha)_p(v_1,\dots,v_k)
       = \alpha_{\phi(p)}\bigl(\phi_{*\,p}(v_1),\dots,\phi_{*\,p}(v_k)\bigr).
     $$

     拉回算子与楔积满足
     $\phi^*(\alpha\wedge\beta)=\phi^*\alpha\wedge\phi^*\beta$。

3. **度量的 pullback（诱导度量）**
   若 $\phi$ 是一个浸入／嵌入，或更一般任何映射，也可对度量做拉回：
   $$
     g' \;=\;\phi^*h
     \;:\;T_pM\times T_pM\;\longrightarrow\;\Bbb R,\quad
     g'(u,v)=h\bigl(\phi_{*\,p}u,\;\phi_{*\,p}v\bigr).
   $$

   $g'$ 是 $M$ 上的一个对称双线性形式；当 $\phi$ 是等距嵌入时，$g'=g$。

#### 4.3.2. Pullback 的代数与微分性质

1. **交换外微分**

   $$
     d\bigl(\phi^*\alpha\bigr)
     =\phi^*(d\alpha),
   $$

   即拉回算子与外微分 $d$ 互换。

2. **与 interior product 对偶**
   若 $X$ 是 $M$ 上向量场，则
   $$
     \iota_X(\phi^*\alpha)
     =\phi^*(\iota_{\phi_*X}\,\alpha).
   $$

3. **复合映射的拉回**
   对 $\psi\colon L\to M$ 也有
   $$
     (\phi\circ\psi)^*
     =\psi^*\circ\phi^*.
   $$

#### 4.3.3. 在黎曼几何中的典型应用

**子流形的诱导度量**：

若 $\iota\colon S\hookrightarrow M$ 是嵌入，则子流形 $S$ 自带的“第一基本形式”正是度量的拉回：

$$
  g_S = \iota^*g.
$$

例如，单位球面 $S^2\subset\Bbb R^3$ 上的度量由 $\iota^*(dx^2+dy^2+dz^2)$ 给出。

**度量测度与积分的变换公式**：

拉回度量后可得到体积形式 $\phi^*(\mathrm{d}V_h)$，并用于积分替换：

$$
\int_M \phi^*(\omega)
  = \int_{\phi(M)} \omega,
$$

若 $\phi$ 是定向同胚，则左边即在 $M$ 上对形式 $\omega$ 拉回后积分。
在局部坐标中，对体积形式的拉回对应雅可比行列式 $\lvert\det d\phi\rvert$。

#### 4.3.4. 例子：圆周映射

令 $\phi\colon S^1\to S^1$ 为“多圈映射”

$$
  \phi(e^{i\theta})=e^{in\theta}.
$$

* **函数拉回**：对 $f(e^{i\theta})=e^{i\theta}$，有 $\phi^*f(e^{i\theta})=e^{in\theta}$。
* **1-形式拉回**：标准角形式 $\omega=d\theta$ 满足
  $\phi^*(d\theta)=d(n\theta)=n\,d\theta$，反映绕圈次数的加权。

#### 4.3.5. 小结

* **Pullback** 提供了一套机械将 $N$ 上的函数、形式、度量等结构“带回”到 $M$ 上的工具；
* 它保留了楔积、外微分、度量内积等代数与微分运算的相容性；
* 在黎曼几何中，最典型的是**诱导度量**与**积分变换公式**，以及子流形几何中第一／第二基本形式的构造基础。

理解 pullback 的这些机制，是在流形间传递微分几何结构，引入子流形、映射几何性质研究（如等距、曲率比较）的前提。



### 4.4. 插入（Contraction，内积在高维流形上的推广）

#### 4.4.1. 坐标环境中的向量场与微分形式

* 局部坐标系是 $(x^1,\dots,x^n)$

* 设微分形式为：

  $$
  \alpha = f(x)\, dx^{i_1} \wedge dx^{i_2} \wedge \cdots \wedge dx^{i_k}
  $$

  也就是说，$\alpha$ 是一个 $k$-形式，它在某些方向上展开出来（比如 $dx \wedge dy$，或 $dx \wedge dz \wedge dt$，等等）。

* 向量场 $X$ 写作：

  $$
  X = a(x)\, \frac{\partial}{\partial x^{j}}
  $$

  ——它是某个方向上的切向量，例如可能是 $\frac{\partial}{\partial x}$，或者某个线性组合。

我们现在要做的，是**把 $X$ 插入 $\alpha$**，变成一个 $(k-1)$-形式。

#### 4.4.2. 插入的几何意义

微分形式的输入是向量。

举个例子：

* $\alpha = dx \wedge dy$，它输入两个向量 $(v_1, v_2)$，输出 $\det(v_1, v_2)$

* 如果你插入一个向量 $X$，就变成了一个 1-形式：

  $$
  (\iota_X \alpha)(v) := \alpha(X, v)
  $$

  即：把第一个槽用 $X$ 填上，剩下一个槽。

所以，插入算子是一个**降阶操作**，几何上对应“在某个方向上切一刀”。

#### 4.4.3. 实例讲解

**示例 1**：

* **1-形式**：设 $\alpha=f\,dx+g\,dy$，向量场 $X=a\,\partial_x+b\,\partial_y$，则
  $$
    \iota_X\alpha
    =f\,a + g\,b
    \;=\;\alpha(X).
  $$

* **2-形式**：令 $\beta=P\,dx\wedge dy$，则

  $$
    \iota_X\beta
    =P\bigl(a\,dy - b\,dx\bigr),
  $$

  因为把 $X=a\,\partial_x+b\,\partial_y$ 插入第一个槽后，剩下的 1-形式是 $P(a\,dy - b\,dx)$。

**示例 2**：

设：

* $\alpha = dx \wedge dy \wedge dz$
* $X = \frac{\partial}{\partial x}$

则插入操作是：

$$
\iota_{\partial/\partial x} (dx \wedge dy \wedge dz) = dy \wedge dz
$$

再来一个：

* $X = \frac{\partial}{\partial y}$，

$$
\iota_{\partial/\partial y} (dx \wedge dy \wedge dz) = - dx \wedge dz
$$

这个符号“负号”的出现，其实就是来源于反对称性 —— 把 $y$ 插进去是“从第二个位置”进来，换到第一位要交换一次，带一个负号。

#### 4. 一般规则（无需任何张量）

现在我们可以总结出**微分形式语言中完全几何化的插入法则**：

> **插入一个向量场 $X$** 到形式 $dx^{i_1} \wedge \cdots \wedge dx^{i_k}$ 中，就是从这个楔积中选一项，把该方向的微分去掉，把向量投影到那一方向的系数乘进去，并加上适当的正负号（看换了几次顺序）。

具体来说，设 $X=\sum_{i} X_i\partial_i$，则
$$
i_X(dx_1\wedge dx_2 \wedge \cdots \wedge dx_n)=\sum (-1)^{i-1} X_i dx_1\wedge \cdots \wedge \hat{dx_i}\wedge \cdots \wedge dx_n
$$
