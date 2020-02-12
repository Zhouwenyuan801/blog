---
title: "ICEPCF"
date: 2020-02-12T16:44:32+08:00
lastmod: 2020-02-12T16:44:32+08:00
draft: false
keywords: ["un","specified"]
description: "There's no description"
tags: ["no-tag"]
categories: ["unarchieved"]
author: "ZZQ"
---

<!--more-->





试试公式功能，$Q_n$，
$$
S_h \times \underbrace{\nabla\cdot \vec{U}}_{convection} = Q
$$


试试代码功能

```c++
#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"   #界面性质
#include "twoPhaseMixtureEThermo.H"  #两相流体热力性质
#include "temperaturePhaseChangeTwoPhaseMixture.H" #关于相变模型
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
volScalarField& T = thermo->T();
    volScalarField& e = thermo->he(); #起别称，方便使用
    e.oldTime(); #create field0Ptr_

    turbulence->validate();#创建后验证湍流模型完整性

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
        rho = alpha1*rho1 + alpha2*rho2;#修正密度项
        runTime.write();
        runTime.printExecutionTime(Info);
    }
    Info<< "End\n" << endl;
    return 0;
}
```

图片，

![奥利给！](https://tva1.sinaimg.cn/large/0082zybply1gbtrgreud3j30bu0fi752.jpg)