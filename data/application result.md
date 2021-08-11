The CpGs found with BY adjust method:

after [2*n/log(n)] SIS between M and X, we get 248 CpGs. After applying SCAD penalized regularization,  we get 25 CpGs. Then we apply Sobel Test to these CpGs, with BY adjust method, we get two CpGs: 

| |cg08530838| cg19757631 |
|:---:|:---:|:---:|
|BY_adj_p|0.04537071 |0.04537071 |
|alpha coef|-0.1352534 |-0.1056452 |
|beta coef|0.07616241 |-0.11213447 |
|alpha\*beta|-0.01030123  |0.01184647 |

```
Bakulski, K. M., Dou, J., Lin, N., London, S. J., & Colacino, J. A. (2019). DNA methylation signature of smoking in lung cancer is enriched for exposure signatures in newborn and adult blood. Scientific reports, 9(1), 1-13.
```

In lung adenocarcinoma samples, there were 14 CpGs associated with current smoking status at a Bonferroni adjusted genome-wide significance level (*P* < 10−7) (Table [3](https://www.nature.com/articles/s41598-019-40963-2#Tab3)), **66 CpGs associated with current smoking stats at FDR significance (Supplementary Table [1](https://www.nature.com/articles/s41598-019-40963-2#MOESM1))**

**Supplementary Table 1.**  DNA methylation sites associated (*FDR* value < 0.05) with smoke exposure (versus never smokers) in The Cancer Genome Atlas lung adenocarcinoma tissue samples. 

Current Smokers

| **CpG** | **Chr** | **Position** | **Annotated Gene** | **Estimated Change in % Meth.** | **Std. Error** | **P-value** | **Average % Meth.** | **Gap Probe** |
| ------- | ------- | ------------ | ------------------ | ------------------------------- | -------------- | ----------- | ------------------- | ------------- |
| cg19757631 | chr1 | 11118889 | SRM  | -12.28 | 2.64 | 4.81E-06 | 60.4 | N    |

```
@article{bakulski2019dna,
  title={DNA methylation signature of smoking in lung cancer is enriched for exposure signatures in newborn and adult blood},
  author={Bakulski, Kelly M and Dou, John and Lin, Nan and London, SJ and Colacino, JA},
  journal={Scientific reports},
  volume={9},
  number={1},
  pages={1--13},
  year={2019},
  publisher={Nature Publishing Group}
}
```

|            | ab_ceof | p value BH |   SE   |  beta   |  alpha  |   low   |   up    |
| :--------- | :-----: | :--------: | :----: | :-----: | :-----: | :-----: | :-----: |
| cg08530838 | -0.0328 |   0.0112   | 0.0099 | 0.2424  | -0.1353 | -0.0521 | -0.0134 |
| cg24720672 | 0.0269  |   0.0151   | 0.0086 | -1.4889 | -0.0181 | 0.0100  | 0.0438  |
| cg19757631 | 0.0296  |   0.0112   | 0.0086 | -0.2806 | -0.1056 | 0.0129  | 0.0464  |
| cg05147638 | 0.0185  |   0.0422   | 0.0070 | 0.4918  | 0.0376  | 0.0047  | 0.0323  |
| cg08636115 | 0.0263  |   0.0152   | 0.0087 | -0.3811 | -0.0690 | 0.0093  | 0.0433  |

