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



## 主文件 interCondensatingEvaporatingFoam.C

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

下面将分别分析。

## 相变模型分析



