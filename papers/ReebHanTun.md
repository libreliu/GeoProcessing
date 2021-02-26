# An Efficient Computation of Handle and Tunnel Loops via Reeb Graphs

[Paper @ ACM](https://dl.acm.org/doi/pdf/10.1145/2461912.2462017)

## 算法

### Handle 和 Tunnel 计算

#### 计算曲面网格 $M$ 的 Reeb Graph，记作 $ Rb_M $

- 计算时采用 $|M| \rightarrow R $ 的一个相对于任意固定平面投影的距离来做高度函数

* 问题：Reeb Graph 如何表示？
  * 可以参考 pxy 提供的文章
#### 计算 $ Rb_M $ 的 "specific cycle basis"

> Rb_M 的 cycle basis 是个啥？
>
> - 就是一个无向图的 cycle basis
> - A fundamental cycle basis may be formed from any spanning tree or spanning forest of the given graph, by selecting the cycles formed by the combination of a path in the tree and a single edge outside the tree. Alternatively, if the edges of the graph have positive weights, the minimum weight cycle basis may be constructed in polynomial time. 
> - https://en.wikipedia.org/wiki/Cycle_basis

对于每个 Rb_M 中的边 e，赋予一个权重 w(e) = h(v_elower)
然后根据这些权重，在 Rb_M 上计算最大权重生成树 T
对每个 Rb_M 中不在 T 上的边诱导一个在 Rb_M 中的 "canonical cycle" c_e

- 具体地说，c_e 是一些 e 唯一确定的环和一些 T 中的边（？）
- 对于任何生成树 T，{c_e} 是 Rb_M 的一个 cycle basis
并且对于这样生成的边权重的 T，有如下性质：
- e 的最低点 p 是 c_e 的最低点，p 是 up-fork saddle

#### 把之前处理的 $ \{c_e\} $ 映射回曲面 M

对于每个 Reeb Graph 的圈 c_e，找一个在曲面 M 上能反映相同路径的一个圈

- 怎么找？
记 \gamma_1 到 \gamma_g 是按最低点高度值升序排列的 g 个圈，鉴于 h 是 morse function，这些圈的最低点高度值应该互不相同，并且每个最低点都是 M 中的 splitting saddle

#### 对于每个圈 $\gamma_i$，我们计算 level set 中的对偶圈 $\bar {\gamma_i}$

- level set 就是那个根据高度产生的等价类中往上偏离一点点的时候的某边的轮廓形成的圈

4. 对每个圈和对偶圈 "perturb“ （往里或者往外”一点点“）
然后分类讨论确定圈是 I 中的还是 O 中的
- 这里似乎只讨论了对偶圈 perturb，非对偶圈怎么搞的（i.e. \gamma_i）？
根据 perturb 后的 configuration 的情况，可以将每种（原来的或者对偶的）圈分成两种：nontrivial in H1(I) 和 nontrivial in H1(O) 的
这样就可以有 \Gamma_1 在 H_1(I) 非平凡，\Gamma_2 在 H_1(O) 非平凡
并且可以证明这两个是 H_1(I) 和 H_1(O) 的 cycle basis
- 但是这个得到的圈不一定是在另一半空间中 trivial 的，为了计算 handle/tunnel loops，我们需要计算 linking number matrix (step 5)
之后在 \Gamma_1 和 \Gamma_2 中分别往里和外面挪动一点（用 level set 的技术），这样的结果和原来不相交，从而为后面算 linking number 做准备

#### 计算 linking number

- 引理：如果某个圈 l 在 \Gamma_1 的每个圈都是 linking number = 0，但是在 \Gamma_2 中至少有一个圈是非 0 的 linking number，那么就是 tunnel loop，handle loop 也类似
在此基础上可以把 L 组成一个矩阵，求出 2g 个满足这种要求的 loop。
这样就可以得到这些 handle 和 tunnel loop。



### 用几何优化的方法优化这组 basis

> 在 $H_1(M)$ 的最短 cycle basis 上的任何圈 $l$ 都是正则 (canonical) 的，即 $ l $ 可以分解为某个 $ e = e(u,v) $ 和两个最短路径 $ \pi(u, w) $ 和 $\pi(v, w)$

