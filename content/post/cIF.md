---

title: "compressibleInterFoam概览"
date: 2020-04-28T13:03:32+08:00
lastmod: 2020-04-28T13:03:32+08:00
draft: false
keywords: ["compressibleInterFoam","OpenFOAM"]
description: "无奈啊"
tags: ["OpenFOAM","code"]
categories: ["OpenFOAM"]
author: "ZZQ"
---

<!--more-->

## 前言

iRMIF本身并不支持可压流体的计算，这点在求解器的描述中已清楚地点明。然而自增压过程天然的是可压流动的范畴，因为气侧的密度、温度会极大的影响封闭体系内的压力。即使如前面博客所述将密度改为温度的函数，也会出现两个问题：1，极易发散。2，压力计算不正常，表现为并无增长或减小，只是在一值附近波动。因此，只能放弃使用iRMIF求解自增压过程。重新看标准求解器，发现有两个重要的选项：1，`compressibleInterFoam`;2，`reactingTowPhaseEulerFoam`。前者是基于`interFoam`发展的可压流动求解器，后者则是另一套思路，对每一相单独求解NS方程，在相间引入拖曳力、升力模型等。虽然Euler方法事实上所需网格明显少于VOF方法，但是因为要引入另外的模化方程，其准确性又变得难以估计。因此，决定在`compressibleInterFoam`的基础上，耦合进相间传热传质模型。

相变模型的耦合同样有两种思路：1. 参考`interCondensatingEvaporatingFoam`的组织架构。2. `icoReactingMultiInterFoam`的组织架构。前面已经分析过，`iCEF`在低温下失效主要是因为其没有使用标准的热物理库接口，同时其处理相变产生的源项时的方法不是很稳定。而`iRMIF`内很多内容是目前不需要的，比如多相多组分间的传质。`cIF`本身使用了热物理类的接口，因此，初定的思路是按照`iCEF`的架构，同时修改一些细节。但是无论如何，要先理清`cIF`求解器的逻辑。



## 求解器主题结构

首先要看的是`cIF`中数据存储结构。在`createFields`中，声明了`twoPhaseMixtureThermo`对象`mixture`（和`interFoam`一致），该对象中储存了几乎所有的两相的物性读取和计算方法。在构建该对象时，根据初始速度、通量调用`phiThermo`，`twoPhaseMixture`，`interfaceProperties`的构造函数，并且构造两个`rhoThermo`对象。之后通过`mixture`对象初始化相分布以及其他的物性。

在主文件中，求解过程也与`interFoam`非常相似，区别是，在`PIMPLE`循环中，调用的是`compressibleAlphaEqnSubCycle.H`以求解相分布。其余温度、压力和速度方程的调用位置均未改变。



## 控制方程

对比`iRMIF`,`UEqn`基本没有区别。`TEqn`中，RHS少了因相变引入的源项。主要区别在于`alphaEqn.H`和`PEqn.H`。下面具体介绍该可压求解器中压力、相分布求解的改变，主要参考[这里](https://www.cfd-online.com/Forums/openfoam-solving/70965-formulation-compressibleinterfoam-2.html)的讨论。

混合物的密度为，
$$
\rho = \alpha_1 \rho1 + \alpha2 \rho2
$$
那么对于每一相的质量连续性方程为，
$$
\frac{d}{dt}(\rho_i \alpha_i )+\nabla\cdot(\rho_i\alpha_iU) = 0
$$
考虑每一个时间步内的密度是常量，有
$$
\frac{d}{dt}(\alpha_i )+\nabla\cdot(\alpha_iU) = -\frac{\alpha_i}{\rho_i}\frac{D}{Dt}\rho_i
$$
考虑每一相的压缩性， $\rho_i= \rho_0 + \psi_i p$，其中$\psi_i$为压缩因子，$p$为压力，上式可以写为式（1），
$$
\frac{d}{dt}(\alpha_i )+\nabla\cdot(\alpha_iU) = -\frac{\alpha_i \psi_i}{\rho_i}\frac{D}{Dt}p
$$
把两相相加，左侧为速度的散度，即式（2）
$$
\nabla\cdot(U)=-(\frac{\alpha_1\psi_1}{\rho_1}-\frac{\alpha_2\psi_2}{\rho_2})
$$


式(1)的RHS可以写为
$$
RHS = \alpha_1 \nabla\cdot(U) + \alpha_1 \dot{g}
$$
其中$\dot{g}$是源项，$\dot{g} = (\psi_2/\rho2 - \psi1/\rho_1)D/Dt(p)$。

因为压力修正过程中用到了连续性方程，所以对应的压力求解也有区别。我们先来看`PEqn.H`。

与`interFoam`或`icoFoam`相同，开始构造了离散的压力泊松方程，

```c++
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi))
    );
    MRF.makeRelative(phiHbyA);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;
```



对应
$$
\bold{HbyA} = \frac{1}{A_p}(-\sum A_N \vec{U}^*_N +\vec{S}^n_P )
$$
具体可以参考这篇[博客](http://dyfluid.com/icoFoam.html)。

之后，区别于不可压求解器，构建了两个部分`p_rghEqnComp1`，`p_rghEqnComp2`分别对应两相的不可压部分。并赋值为，

```c++
        p_rghEqnComp1 =
            pos(alpha1)
           *(
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f)
                  - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
                )/rho1
              - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh))
            );

        p_rghEqnComp2 =
            pos(alpha2)
           *(
               (
                   fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f)
                 - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
               )/rho2
             - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
             + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh))
            );
```

定义中括号里的对应连续性方程，收敛时应为0。后面的`correction()`表示为新旧两个值得差，新值对应的是式（1），收敛时也为0。因此最后收敛时只剩下$p_{old}$这一项数学表达为
$$
p_{\rho gh,com1} =\frac{1}{\rho_1}(\frac{d}{dt}(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1U))-\frac{d}{dt}(\alpha_1)-\nabla\cdot(\alpha_1U)+\frac{\alpha_1\psi_1}{\rho_1}\frac{D}{Dt}(p-p_{old})
$$

$$
p_{\rho gh,com1} =-\frac{\alpha_1\psi_1}{\rho_1}\frac{D}{Dt}p_{old}
$$

这其实就是式1的RHS。

随后，两部分可压部分加上不可压部分后共同求解，关于压力修正方程的细节可以参考这篇[博客](https://www.cfd-online.com/Forums/openfoam-solving/196243-understanding-terms-compressible-peqn.html)，

```c++
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );
```

在最后一步迭代的时候，赋值给`dgdt`

```c++
            dgdt =
            (
                alpha1*(p_rghEqnComp2 & p_rgh)
              - alpha2*(p_rghEqnComp1 & p_rgh)
            );
```

这里网上给出的结论都是对应的，
$$
\dot{g} = \alpha_1\alpha_2(\frac{\psi_2}{\rho_2}-\frac{\psi_1}{\rho_1})\frac{Dp}{Dt}
$$
这里`p_rghEqnComp1 & p_rgh`起到的作用怎么也查不到，只能先默认是对的。

到这里压力的部分就结束了。

随后，在`compressibleAlphaEqnSubcycle.H`中，调用了`alphaEqn.H`，而`alphaEqn.H`又调用了`alphaSuSp.H`。在`alphaSuSp.H`中，

```c++
forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
}
```

即将源项按对角占优原则处理，在相方程中，求解的是，

```c++
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
       // - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
       //           + fvc::div(phiCN), alpha1)
         ==
            Su + fvm::Sp(Sp + divU, alpha1)
        );
```

简单推导后，代入式（2），即可得到该式即对应式（1）。至此，控制方程就分析完了。

如果需要在该求解器中耦合传热传质，要考虑的事情有：

1. 在mixture类中植入相间传质速率的计算方法，可以从`interfaceComposition`着手，也可以简单地使用`iCEF`中的方法

2. 类似的，参考`iRMIF`或`iCEF`，做出相间传热的计算方法。

3. 在`TEqn`，`UEqn`，和`alphaSuSp`中，将相变引入的源项做相应的处理。

   对应方程的修改应该追溯到
   $$
   \frac{d}{dt}(\rho_1\alpha_1) + \nabla\cdot(\rho_1\alpha_1U) = \dot{m}
   $$
   