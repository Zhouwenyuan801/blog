---
title: "interCondensatingEvaporatingFoam笔记"
date: 2020-02-12T16:44:32+08:00
lastmod: 2020-02-12T16:44:32+08:00
draft: false
keywords: ["E&C","OpenFOAM"]
description: "OpenFOAM中蒸发冷凝求解器笔记"
tags: ["E&C","OpenFOAM","code"]
categories: ["OpenFOAM"]
author: "ZZQ"
typora-copy-images-to: ipic
---

<!--more-->

## 总览

`interCondensatingEvaporatingFoam`是自ESI-OpenCFD发布的`OpenFOAM-1806`版本以来新增的求解器，其描述为：

>​    Solver for two incompressible, non-isothermal immiscible fluids with
>
>​    phase-change (evaporation-condensation) between a fluid and its vapour.
>
>​    Uses a VOF (volume of fluid) phase-fraction based interface capturing
>
>​    approach.
>
>
>
>​    The momentum, energy and other fluid properties are of the "mixture" and a
>
>​    single momentum equation is solved.
>
>
>
>​    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

其中，相方程的建立与求解与`interFoam`相同，涉及相变的内容多沿袭自`interPhaseChangeFoam`。对于蒸发冷凝过程，使用`temperaturePhaseChangeTwoPhaseMixture`类描述。在该求解器在早期版本（v1806）中，流体表面张力、温度与常温流体相差较大时计算较易发散。在最近的版本（v1912）中，模型中的表面张力模型、能量方程进一步更新，修复了前述问题。因此，该求解器目前已较为完善。遗憾的事，在OpenFOAM Fundation发行的`OpenFOAM-6`及后续版本中并没有该求解器。好在学校的超算预置的版本是v1812。

后续将记录分析该求解器的过程，很多内容都是较基础的知识，这部分知识只作简单总结。



## 主文件

该求解器引用的头文件如下，除必要的引用外，还有界面属性类、热物性类和相变混合物体系类的引用。

```c++
#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"   #界面性质类
#include "twoPhaseMixtureEThermo.H"  #热物性类
#include "temperaturePhaseChangeTwoPhaseMixture.H" #相变混合物体系类
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
```

主函数为

```c++
{
  #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"  #初始化场
    #include "createFvOptions.H"
    #include "createTimeControls.H" #建立时间序列
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volScalarField& T = thermo->T();  #获得温度场引用
    volScalarField& e = thermo->he(); #起别称，方便使用
    e.oldTime(); #create field0Ptr_建立时间序列

    turbulence->validate(); #创建后验证湍流模型完整性

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "createTimeControls.H"
        #include "CourantNo.H" #计算Co数
        #include "setDeltaT.H" #根据Co修正时间步长
        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop()) #PIMPLE循环
        {
            #include "alphaControls.H"
            surfaceScalarField rhoPhi  #建立质量通量
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, Zero)
            );
            mixture->correct();#修正混合流体物性
            #include "alphaEqnSubCycle.H"
            solve(fvm::ddt(rho) + fvc::div(rhoPhi)); #体积分数advection
            #include "UEqn.H"  #求解速度
            #include "eEqn.H"  #根据速度修正能量
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"  #压力修正
            }
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        rho = alpha1*rho1 + alpha2*rho2; #修正密度项
        runTime.write();
        runTime.printExecutionTime(Info);
    }
    Info<< "End\n" << endl;
    return 0;
}
```

该求解器构建中比较重要的是`temperaturePhaseChangeTwoPhaseMixture`类，其负责处理相变量的计算。该类继承自`IOdictionary`，包含两个成员对象，`twoPhaseMixtureEThermo`和`mesh`。

![image-20200212183756118](https://tva1.sinaimg.cn/large/0082zybply1gbtswsmk8tj30x20i4did.jpg)

`twoPhaseMixtureEThermo`主要包含两相流体的热物性质，`mesh`包含流场中的网格信息。`temperaturePhaseChangeTwoPhaseMixture`则将热物理性质映射至网格上。`twoPhaseMixtureEThermo`继承自`basicThermo`和`thermoIncompressibleTwoPhaseMixture`，前者是标准类，包含焓、温度、压力等成员变量，后者包含潜热、比热容、热导率等与传热相关的物性。

在该求解器中，重要的内容有：

- 混合流体物性计算和更新
- 蒸发冷凝模型的耦合
- 表面张力模型

下面将分别分析上述内容。

## 相变模型分析

本求解器的近亲`interPhaseChangeFoam`处理的是空化过程中的相变，这里处理的是由热量传递引起的相变。

相较空化，除同样关注压力外，界面附近的温度也是模型的重要参量。本求解器的相变过程由`temperaturePhaseChangeTwoPhaseMixture`类处理。先分析该类的头文件。

```c++
    //- Runtime type information,作为父类
    TypeName("temperaturePhaseChangeTwoPhaseMixture");

    // Declare run-time constructor selection table，构建RTS表，key值为模型名，value为指针
        declareRunTimeSelectionTable
        (
            autoPtr,
            temperaturePhaseChangeTwoPhaseMixture,
            components,
            (
                const thermoIncompressibleTwoPhaseMixture& mixture,
                const fvMesh& mesh
            ),
            (mixture, mesh)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model，New函数接受两个参数
        static autoPtr<temperaturePhaseChangeTwoPhaseMixture> New
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );
```

选择器的实现为，

```c++
Foam::autoPtr<Foam::temperaturePhaseChangeTwoPhaseMixture>  #返回指向相变模型的指针
Foam::temperaturePhaseChangeTwoPhaseMixture::New
(
    const thermoIncompressibleTwoPhaseMixture& thermo,
    const fvMesh& mesh
)
{
    IOdictionary phaseChangePropertiesDict    #读取constant目录下的phaseChangeProperties
    (
        IOobject
        (
            "phaseChangeProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType   #读取相变模型字符串
    (
        phaseChangePropertiesDict.get<word>("phaseChangeTwoPhaseModel")
    );

    Info<< "Selecting phaseChange model " << modelType << endl;

    auto cstrIter = componentsConstructorTablePtr_->cfind(modelType); #在表中寻找对应key值

    if (!cstrIter.found()) #若没有输出表中所有key值
    {
        FatalErrorInFunction
            << "Unknown temperaturePhaseChangeTwoPhaseMixture type "
            << modelType << nl << nl
            << "Valid temperaturePhaseChangeTwoPhaseMixture types :" << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<temperaturePhaseChangeTwoPhaseMixture>
        (cstrIter()(thermo, mesh));
}
```

关于RTS，有几篇博客写的很好[[1]](http://xiaopingqiu.github.io/2016/03/12/RTS1/)[[2]](https://marinecfd.xyz/post/openfoam-runtime-selection/)，本有一篇英文文章更详细些，可惜站点挂掉了。

该类内关于相变的函数如下，英文本身的注释很不清晰，这里按我的理解重新解释，

```c++
//相方程中质量传递源项，纯虚
virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;
//相变模型直接计算得到的相变速率，kg/m^3，在压力方程中使用，纯虚
virtual Pair<tmp<volScalarField>> mDot() const = 0;
//没弄懂
virtual Pair<tmp<volScalarField>> mDotDeltaT() const = 0;
//对mDotAlphal包装
Pair<tmp<volScalarField>> vDotAlphal() const;
//对mDot包装
Pair<tmp<volScalarField>> vDot() const;
//更新相变场
virtual void correct() = 0;
//读取相关参数
virtual bool read();
```



理解该模型需要从PISO算法和相方程开始：

在PISO算法中，压力修正方程右侧原本为0，但因为质量守恒方程有源项，所以方程更新为，
$$
\sum_{f\in\partial \Omega_i}\left ( \frac{1}{A_P}\right )\left ( \nabla^\perp_f p^{m+1}_d\right ) |S_f| - \sum_{f\in\partial \Omega_i} \phi^r_f = \frac{\dot{m}}{\rho_l}- \frac{\dot{m}}{\rho_v}
$$
而考虑到相变的相输运方程为，
$$
\frac{\partial \alpha_l}{\partial t} + \nabla \cdot (\alpha_l \mathbf U) =\alpha_l\nabla\cdot\mathbf U-\alpha_l \left ( \frac{\dot{m}}{\rho_l}- \frac{\dot{m}}{\rho_v} \right ) + \frac{\dot{m}}{\rho_l}
$$
右侧前两项在不可压流体中为0，计算时增加可以提高数值稳定性。MULES对该方程有进一步的修正，并没有直接体现在这个方程内。

以派生类`constant`为例介绍该模型的求解过程。`constant`类实际即为`Lee`模型，该模型中相变量表示为，
$$
\dot{m}_l=r_e\alpha_l\rho_l\frac{T-T_\mathrm{sat}}{T_\mathrm{sat}}
$$
`constant`中`mDot()`函数定义为，

```c++
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDot() const
{

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*limitedAlpha2*max(TSat - T.oldTime(), T0),
        coeffE_*mixture_.rho1()*limitedAlpha1*max(T.oldTime() - TSat, T0)
    );
}
```

这里将饱和温度同时考虑进系数中，关于该类模型变种太多，这里不赘述。与`mDot()`类似的`mDotAlphal()`函数内容为，

```c++
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotAlphal() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*max(TSat - T.oldTime(), T0),
       -coeffE_*mixture_.rho1()*max(T.oldTime() - TSat, T0)
    );
}
```

该函数相较`mDot()`仅少了对应相的质量分数，并且蒸发项多了负号。初看十分不理解，要在具体方程中还原才能看清。在`alphaEqn.H`和`pEqn.H`中，相方程和压力修正方程定义为，

```c++
Pair<tmp<volScalarField>> vDotAlphal =
        mixture->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
        );
```

```c++
    Pair<tmp<volScalarField>> vDot = mixture->vDot();
    const volScalarField& vDotc = vDot[0]();
    const volScalarField& vDotv = vDot[1]();

        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          - (vDotc - vDotv)
        );
```

回看`vDot()`和`vDotAlphal()`的定义，

```c++
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixture::vDot() const
{
    dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());
    Pair<tmp<volScalarField>> mDot = this->mDot();

    return Pair<tmp<volScalarField>>(pCoeff*mDot[0], pCoeff*mDot[1]);
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff
    (
        1.0/mixture_.rho1() - mixture_.alpha1()
       *(1.0/mixture_.rho1() - 1.0/mixture_.rho2())
    );

    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}
```

较易理解的事压力方程，因为$\dot{m}=\mathrm{mDot[0]-mDot[1]}$，`pEqn`即为方程（1）

而带入表达式，`alphalEqn`右侧为
$$
\mathrm{fvm::Sp(vDotvmcAlphal, alpha1) + vDotcAlphal} 
$$

$$
[\frac{1}{\rho_l} - \alpha_l (\frac{1}{\rho_l}-\frac{1}{\rho_g})]\cdot (S_g -S_l )\cdot \alpha_l +[\frac{1}{\rho_l} - \alpha_l (\frac{1}{\rho_l}-\frac{1}{\rho_g})]\cdot S_l
$$

$$
[\frac{1}{\rho_l} - \alpha_l (\frac{1}{\rho_l}-\frac{1}{\rho_g})] \cdot [(S_g -S_l)\alpha_l +S_l]
$$

$$
[\frac{1}{\rho_l} - \alpha_l (\frac{1}{\rho_l}-\frac{1}{\rho_g})] \cdot [S_g \alpha_l +S_l \alpha_g]
$$

考虑$S_g$函数表达式，可以看出`alphalEqn`即为相输运方程。

## 表面张力模型

前述中存在的问题**有可能**与表面张力模型的更新有关。表面张力影响速度方程，在`UEqn`中，速度方程定义为，

```c++
fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
```

表面张力项为`interface.surfaceTensionForce()`。对应的动量方程为，
$$
\frac{\partial\rho \mathbf U}{\partial t} + \nabla\cdot(\rho\mathbf U \mathbf U)-\nabla\cdot\tau=-\nabla p+\rho\cdot h\nabla\rho +F_\mathrm{st}
$$
该求解器使用CSF模型，将表面张力表示为[[3]](http://dyfluid.com/interFoam.html)，
$$
F_\mathrm{st} = \sigma \kappa \nabla\alpha_l
$$
其中，$\sigma$是表面张力系数， $\kappa$ 是界面曲率，
$$
\kappa = -\nabla \cdot \mathbf n
$$

$$
\mathbf n = \frac{\nabla\alpha_l }{|\nabla \alpha_l |}
$$

$\mathbf n$ 是自由界面的法相。表面张力计算、更新由`interfaceProperties`类实现，由速度、相分布场构建，其成员变量有

```c++
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension, \sigma,这里使用了RTS
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal, n
        const dimensionedScalar deltaN_;
				//- Reference to the field
        const volScalarField& alpha1_;
        const volVectorField& U_;
				// \vec{n}_f
        surfaceScalarField nHatf_;
				// curvation, \kappa
        volScalarField K_;
```

主要的成员函数有，

```c++
private:
        //- Re-calculate the interface curvature
        void calculateK();
public:
        tmp<volScalarField> sigmaK() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;
        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void correct();
```

先看该类的构造函数，

```c++
Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(alpha1.mesh().V()))
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    )
{
    calculateK();
}
```

考虑`creatFields.H`中本例的构建，

```c++
Info<< "Creating temperaturePhaseChangeTwoPhaseMixture\n" << endl;
autoPtr<temperaturePhaseChangeTwoPhaseMixture> mixture =
    temperaturePhaseChangeTwoPhaseMixture::New(thermo(), mesh);

interfaceProperties interface
(
    alpha1,
    U,
    thermo->transportPropertiesDict()
);
```

初始构建时使用液相体积分数，速度场，以及混合物类中的输运字典，并调用`calculateK()`初始化曲率。表面张力系数由New出来的指针返回，可以是定值，也可以是方程拟合。使用方程的派生类目前还未实施，仅能通过编译。

```c++
void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    //离散面法向量
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
  
    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
		// \nabla n_f / | \nabla n_f |
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);
}
```

这里有一个容易误解的地方，`nHatfv`即是面上插值的自由界面法向，得到后再进行`fvm::div`才得到曲率。这里的`div`不是求散度，而是求和[[4]](http://xiaopingqiu.github.io/2015/05/17/OpenFOAMcode1/) [[5]]([http://www.cfd-china.com/topic/1474/openfoam%E5%AF%B9%E6%A0%87%E9%87%8F%E6%B1%82%E6%95%A3%E5%BA%A6/3](http://www.cfd-china.com/topic/1474/openfoam对标量求散度/3))，
$$
\sum_f \mathbf n_f \cdot S_f = \oint_{\partial V} \mathbf n \cdot dS = \int_V \mathbf n dV
$$
计算出曲率后，表面张力有`surfaceTensionForce()`计算，

```c++
Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}
```

求解中每一步后调用`correct()`对表面张力进行修正，其实就是重新计算界面曲率。

```c++
void Foam::interfaceProperties::correct()
{
    calculateK();
}
```

除此之外，考虑边界条件对表面张力的修正及接触角计算等这里不再赘述。



## 混合物物性

这一部分写的比较乱，官方很多函数还没有写完。密度、比热容、热导率、热扩散系数和内能均有相系数加权求和得到。能量方程中求解的是内能，得到内能场后反推温度场，温度与内能的关系为，
$$
e = (\alpha_1\rho_1e_1 + \alpha_2\rho_2e_2) / (\alpha_1\rho_1 + \alpha_2\rho_2)
$$
其中，
$$
e_1 = Cv_1(T - T_\mathrm{sat}) + Hv_1
$$

$$
e_2 = Cv_2(T - T_\mathrm{sat}) + Hv_2
$$



因此，温度场计算为
$$
T =
        \frac{
            e (\alpha_1\rho_1 + \alpha_2\rho_2)
         -  (\alpha_1\rho_1Hv_1 + \alpha_2\rho_2Hv_2)
       }
       {\alpha_1\rho_1*Cv_1 + \alpha_2\rho_2Cv_2}
       + T_\mathrm{sat};
$$


对应在`twoPhaseMixtureEThremo::correct()`中，定义为：

```c++
    T_ =
        (
            (e_*(alpha1Rho1 + alpha2Rho2))
         -  (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
        )
       /(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
       + TSat_;

    T().correctBoundaryConditions();
```



至此，该求解器中重要的部分就大体分析完了。