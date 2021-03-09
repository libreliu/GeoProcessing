# Greedy optimal homotopy and homology generators

这篇论文回答了

How does one compute the shortest system of loops, among all systems, relaxing the homotopy condition?

Colin and Verdiere 他们做的是在一个同伦群中的 shortest system of loops

这篇文章是以某个 basepoint 开始的 2g 个 loop，不一定同伦等价

类比 Steiner tree 和 MST
Steiner tree
- https://en.wikipedia.org/wiki/Steiner_tree_problem
- 相当于可以自己“加点”的MST

special case of one-vertex cut graphs, which they called systems of loops

cut locus 是 M 中点的集合的闭包，其中至少有两个从基准点 x 的最短路径

https://en.wikipedia.org/wiki/Cut_locus


对于光滑表面，cut locus 一个有限图到 M 的嵌入，并且 M 是 M \ {x} 的 deformation retract

matroid: 集合 X 的非空子集集合 "independent set"
1. 每个独立集的子集都是独立集
2. 如果独立集 A 和 B 满足 |A| > |B|, 那么存在 a \in A \ B 使得 B \union \{a\} 也是独立集

最大独立集被称为 matroid 的基
每个 matroid 的基都有相同的 cardinality

寻找最短 homology basis 是一个 matroid 优化问题
但是，寻找最优 homotopy basis 不能简单套用在 matriod 优化框架中

定义 partial homotopy basis 为生成基本群的子群 G 的 L 的最小有限子集，但是 2-manifold 的 partial homotopy basis 不一定是 matroid

> greedy homotopy basis：
>
> 对于每个 $ i $，$ \gamma_i $ 是使得 $ M \ (\gamma_1 \cup ... \cup \gamma_{i-1} \cup l) $ 连通的最短圈 $ l $ 



reduced cut locus $ \Phi $, cut path $ \phi \in \Phi $ 

$ \sigma(c, \phi) $ 是含点 $ c $ （$ c \in \phi $）的最短的不可收缩圈，定向为圈从左向右穿过 cut path $ \phi $ 

$ \sigma(b, \phi) $ 是含点 $ b $ （$ b $ 靠近 $\phi $）的最短的不可收缩圈，且与某 $ \sigma(c, \phi) $ 同伦等价

$ \sigma(\phi) $ 是对 $ \phi $ 的闭包中的所有点 $ c $ 形成的 $ \sigma(c, \phi) $ 中最短的那个圈

> Q: 这里的闭包是什么？对什么运算封闭？

引理 3.3 每个 greedy homotopy basis 中的圈都对于某些 cut path $ \phi $ 有形式 $ \sigma(\phi) $



combinatorial surface: path 只能在 1-skeleton 上

piecewise-linear surface: path 可以在面内部



## 对于给定 Combinatorial Surface 计算 Greedy Homology Basis

1. 对于任意顶点 $ x $，计算 reduced cut locus $ \Phi(x) $ 和 greedy homotopy basis $ \gamma_1(x), ..., \gamma_{2g}(x) $
   1. 