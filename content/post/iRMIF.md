---

title: "IRMIF中蒸发冷凝模型实现方式"
date: 2020-02-27T13:03:32+08:00
lastmod: 2020-02-27T13:03:32+08:00
draft: false
keywords: ["E&C","OpenFOAM"]
description: "icoReactingMultiphaseInterFoam中蒸发冷凝模型学习记录"
tags: ["E&C","OpenFOAM","code"]
categories: ["OpenFOAM"]
author: "ZZQ"
---

<!--more-->

## 总览

`icoReactingMultiphaseInterFoam`是相较iCEF更复杂也更全面的一个求解器。除两相间的换热外，该求解器还可以计算多相间的传热传质（比如气相中氦气的增压过程）。相对应的，求解器涉及到的功能更多，主要是做化学的人在维护，但很多功能现在还用不到，所以这次并没有选用。但是这个求解器包含了实际上是`Schrage`模型的相变模型。非常可以参考这种规范性的写法。下面记录分析思路。

## 求解器整体框架分析

先看主文件，主文件中特殊的地方有：

- 似乎所有跟混合物、界面有关的求解都通过`fluid`完成，例如在每个PISO循环前有`fluid.solve()     `,`rho=fluid.rho()`更新混合体系状态
- 除温度、压力和速度方程外，求解器另外求解了混合物浓度方程`YEqn`。该方程中也包含了相方程中源项的计算：

```c++
// Calculate mass exchange for species consistent with
// alpha's source terms.
fluid.massSpeciesTransfer(phase, Sus[i], Sps[i], Y[i].name());
```

在压力方程中，和`iCEF`相同，也有相变的计算，

```c++
tmp<volScalarField> tdmdt12(fluid.dmdt(key12))

  p_rghEqn +=
dmdtNet*
(
- fluid.coeffs(phase1.name())
+ fluid.coeffs(phase2.name())
);
```

其中，`key12`是两相体系组成的键值对，而`fluid.dmdt()`是相变模型返回的相变速率。

那么问题的关键就是弄懂`fluid`到底是什么东西，以及相变模型是如何耦合的。

在`createFileds.H`中发现，`fluid`是指向`multiphaseSystem`类的指针指向对象的引用。

```c++
Info<< "Creating multiphaseSystem\n" << endl;
autoPtr<multiphaseSystem> fluidPtr = multiphaseSystem::New(mesh);
multiphaseSystem& fluid = fluidPtr();
```

这种调用方式说明这其实是RTS动态选取的类，在后文中有介绍。`multiphaseSystem.h`文件在`phaseSystem`文件夹内，而`multiphaseSystem`也是前者的子类。下面先看`phaseSystem`类的构造，其继承自`basicThermo`和`compressibleTransportModel`，并包含有两相体系表、相模型表、表面张力模型表、界面多孔模型表等Hash表。这些表记录了特定体系中使用的具体方法。

## 模型相变框架

`multiphaseSystem`类中多定义了

```c++
//- Maximum volumen rate change
dimensionedScalar ddtAlphaMax_;
//- Compression fluxed for phases
compressionFluxTable limitedPhiAlphas_;
//- Su phase source terms
SuSpTable Su_;
//- Sp phase source terms
SuSpTable Sp_;

void calculateSuSp();
virtual void solve();
```

其中，`Su`，`Sp`为相方程源项，不同于`iCEF`，本求解器将相方程改写为，
$$
\frac{\partial \alpha_l}{\partial t} + \nabla \cdot (\alpha_l \mathbf U) = S_u + \alpha_l S_p
$$


其中，
$$
S_u = \frac{\dot{m}}{\rho_l}
$$

$$
S_p = \frac{\dot{m}}{\rho_l-\rho_v}
$$

同时，该类还增加了主文件中出现的`solve()`函数和看起来像是计算源项的`calculateSuSp()`函数。

`solve()`函数定义了相方程，并完成了MULES修正。

```c++
            fvScalarMatrix alpha1Eqn
            (
                fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phi,
                    upwind<scalar>(mesh, phi)
                ).fvmDiv(phi, alpha1)
              ==
                 Su + fvm::Sp(Sp, alpha1)
            );

            alpha1Eqn.solve();
```

可见，这里将`Sp`项隐式处理进系数矩阵中。而`iCEF`中是以显式源项的形式出现，目前不知道这两种处理方式的差异有多大，但联想到`iCEF`中的骚操作，个人更倾向使用这种方法。在求解相方程前，该函数也调用了`calculateSuSp()`方法。所以该函数完成的任务是求解每一个两相体系中的源项，并完成相迁移的计算。



接下来看`calculateSuSp()`方法。因为涉及到多相体系，这里有很多遍历语句和判断，这里略去不表，只分析其中的一种情况。重要的是更新源项的内容。

```c++
        tmp<volScalarField> tdmdt12(this->dmdt(key12));
        const volScalarField& dmdt12 = tdmdt12();
        tmp<volScalarField> tdmdt21(this->dmdt(key21));
        const volScalarField& dmdt21 = tdmdt21();
        volScalarField::Internal& SpPhase1 = Sp_[phase1.name()];

        volScalarField::Internal& SuPhase1 = Su_[phase1.name()];

        volScalarField::Internal& SpPhase2 = Sp_[phase2.name()];

        volScalarField::Internal& SuPhase2 = Su_[phase2.name()];

        const volScalarField dmdtNet(dmdt21 - dmdt12);

        const volScalarField coeffs12(coeffs1 - coeffs2);
        forAll(dmdtNet, celli)
        {
            scalar dmdt21 = dmdtNet[celli];
            scalar coeffs12Cell = coeffs12[celli];

            scalar alpha1Limited = max(min(alpha1[celli], 1.0), 0.0);

            // exp.
            SuPhase1[celli] += coeffs1[celli]*dmdt21;

            if (dmdt21 > 0)
            {
                if (coeffs12Cell > 0)
                {
                    // imp
                    SpPhase1[celli] -= dmdt21*coeffs12Cell;
                }
                else if (coeffs12Cell < 0)
                {
                    // exp
                    SuPhase1[celli] -=
                        dmdt21*coeffs12Cell*alpha1Limited;
                }
              ...
              ...
```

可见，这里是先获得了`dmdtNet`即$\dot{m}$ ，随后判断相变的方向再按公式(2-3)更新源项。那么相变模型的选取将最终影响到这个`dmdt()`这个函数，但这个函数在`multiphaseSystem.C`中并没有定义，联想到前面提到该类其实是动态选取的，在`multiphaseSystems.C`中发现，

```c++
    typedef
        MassTransferPhaseSystem<multiphaseSystem> massTransferMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        massTransferMultiphaseSystem,
        dictionary,
        massTransferMultiphaseSystem
    );
```

可以看到，这里将`MassTransferPhaseSystem<multiphaseSystem>`加入到了类的RTS表中。这个类是以`multiphaseSystem`为模版参数的派生类，其内容包括

```c++
        //- Overall inter-phase mass transfer rates [Kg/s]
        dmdtTable dmdt_;

        //- Mass transfer models
        massTransferModelTable massTransferModels_;
        //- Return total interfacial mass flow rate
        tmp<volScalarField> dmdt(const phasePairKey& key) const;

        //- Return the heat transfer matrix and fill dmdt for phases
        virtual  tmp<fvScalarMatrix> heatTransfer(const volScalarField& T);

        //- Calculate mass transfer for species
        virtual void massSpeciesTransfer
        (
            const phaseModel& phase,
            volScalarField::Internal& Su,
            volScalarField::Internal& Sp,
            const word speciesName
        );
```

两个表对应各个相对的传质速率和选用的相变模型，`massSpeciesTransfer`如前所述返回的是相方程的RHS，`dmdt()`只是数据接口，真正起作用的是`heatTransfer()`方法，该方法在`TEqn.H`中被调用。那么这个函数起到的作用应该是返回各个计算单元内与相变有关的传热量，同时计算更新传质速率。这应该是本求解器关于相变的重点。因为代码比较长，所以大部分标注在代码上。

```c++
template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    const volScalarField& T
)
{
  //构建增加在TEqn中的Matrix类
    tmp<fvScalarMatrix> tEqnPtr
    (
        new fvScalarMatrix(T, dimEnergy/dimTime)
    );

    fvScalarMatrix& eqn = tEqnPtr.ref();
	//循环所有的两相体系
    forAllConstIters(this->phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != this->phaseModels_.end(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                // Net mass transfer from k to i phase
                tmp<volScalarField> tdmdtNetki
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "tdmdtYki",
                            this->mesh().time().timeName(),
                            this->mesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, Zero)
                    )
                );
                volScalarField& dmdtNetki = tdmdtNetki.ref();

							//如果目前两相体系有定义相变模型的话
                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    // Explicit temperature mass transfer rate, Kexp()
                    tmp<volScalarField> Kexp =
                        interfacePtr->Kexp(interfaceCompositionModel::T, T);

                    if (Kexp.valid())
                    {
                      //更新类私有数据dmdt_
                        dmdtNetki -= Kexp.ref();
                        *dmdt_[keyik] = Kexp.ref();
                    }
                }

                // Looking for mass transfer in the other direction (k to i)
                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    // Explicit temperature mass transfer rate
                    const tmp<volScalarField> Kexp =
                        interfacePtr->Kexp(interfaceCompositionModel::T, T);

                    if (Kexp.valid())
                    {
                        dmdtNetki += Kexp.ref();
                        *dmdt_[keyki] = Kexp.ref();
                    }

                }

                word keyikName(phasei.name() + phasek.name());
                word keykiName(phasek.name() + phasei.name());

                eqn -=
                  //更新传热方程，每一个两相体系都有一项
                    (
                        dmdtNetki
                        *(
                            calculateL(dmdtNetki, keyik, keyki, T)
                            - (phasek.Cp() - phasei.Cp())
                            * constant::standard::Tstd
                        )
                    );
            }
        }
    }
    return tEqnPtr;
}
```

calculateL()的作用只是返回潜热，对于热量的计算与`iCEF`相似。可见，这里最终表示传质速率的是`interfaceComposition`类的`Kexp()`方法。`massTransferModelTable`存储的也正是这个类的对象。

##相变模型选取

注意，源代码在`massTransferModels`文件夹内有`interfaceCompositionModel`和`InterfaceCompositionModle`，拷贝到不区分大小写的系统会损失部分代码。这里有一些细节同样不考虑，只看`Kexp()`。作为模版类，`interfaceCompositionModel`只提供纯虚的`Kexp()`方法，派生出`Lee`，`kineticGassEvaporation`（其实是简化后的Schrage）两种相变模型。下面只看Schrage，该模型表示相变速率为：
$$
\dot{m} = \frac{2}{2-\sigma_c} \sqrt{\frac{M}{2\pi R}} \left[\sigma_c \frac{p_g}{\sqrt{T_{g,sat}}} - \sigma_e \frac{p_l}{\sqrt{T_{l,sat}}} \right ]
$$
这里使用的公式是：
$$
\dot{m} = \frac{2}{2-\sigma_c} \sqrt{\frac{M}{2\pi R T_{sat}^3}} h_{fg} (T-T_{sat})\frac{\rho_l\rho_g}{\rho_l - \rho_g}
$$
依据是：Tanasawa, Advances in condensation heat transfer, in: J.P. Hartnett, T.F. Irvine (Eds.), Advances in Heat Transfer, Academic Press, San Diego, 1991. 

这是一种简化版的Schrage。

```c++
template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field)
{
    if (this->modelVariable_ == variable)
    {
      //这里获得的实际上是相分数场
        const volScalarField& to = this->pair().to();

        const volScalarField& from = this->pair().from();

        const fvMesh& mesh = this->mesh_;
		//上一步的温度场
        const volScalarField& T =
            mesh.lookupObject<volScalarField>("T").oldTime();
		//HK方程前面的系数 \sqrt{\frac{M}{2\pi R T^3}
        const dimensionedScalar HerztKnudsConst
        (
            sqrt
            (
                Mv_
               /2.0
               /constant::physicoChemical::R
               /mathematical::pi
               /pow3(Tactivate_)
            )
        );

        word fullSpeciesName = this->transferSpecie();
        auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));

        tmp<volScalarField> L = this->L(speciesName, field);
				// \nabla \alpha_l
        const volVectorField gradFrom(fvc::grad(from));
        const volVectorField gradTo(fvc::grad(to));
				// 界面密度，可以用来计算界面面积
        const volScalarField areaDensity("areaDensity", mag(gradFrom));
			
        const volScalarField gradAlphaf(gradFrom & gradTo);
			
        volScalarField Tmask("Tmask", from*0.0);
			//Mask场，在界面位置为1
        forAll(Tmask, celli)
        {
            if (gradAlphaf[celli] < 0)
            {
                if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
                {
                    {
                        scalar alphaRes = 1.0 - from[celli] - to[celli];
                        if (alphaRes < alphaRestMax_)
                        {
                            Tmask[celli] = 1.0;
                        }
                    }
                }
            }
        }
			//Schrage中的系数
        tmp<volScalarField> tRhom
        (
            new volScalarField
            (
                IOobject
                (
                    "trhom",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimDensity, Zero)
            )
        );
        volScalarField& rhom = tRhom.ref();
      
        tmp<volScalarField> tTdelta
        (
            new volScalarField
            (
                IOobject
                (
                    "trhom",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimTemperature, Zero)
            )
        );
        volScalarField& tDelta = tTdelta.ref();
			//将相变区域限制在界面附近
        if (sign(C_.value()) > 0)
        {
            rhom =
                this->pair().to().rho()*this->pair().from().rho()
              / (this->pair().from().rho() - this->pair().to().rho());

            tDelta = max
            (
                (T*Tmask - Tactivate_),
                dimensionedScalar("T0", dimTemperature, Zero)
            );
        }
        else
        {
            rhom =
                this->pair().to().rho()*this->pair().from().rho()
              / (this->pair().to().rho() - this->pair().from().rho());

            tDelta = max
            (
                Tmask*(Tactivate_ - T),
                dimensionedScalar("T0", dimTemperature, Zero)
            );
        }
			//Schrage模型给出的相变速率，单位为kg/m^2 s
        volScalarField massFluxEvap
        (
            "massFluxEvap",
            2*mag(C_)/(2 - mag(C_))
          * HerztKnudsConst
          * L()
          * rhom
          * tDelta
        );

        // 'from' phase normalization
        // WIP: Normalization could be convinient for cases where the area were
        // the source term is calculated is uniform
      //这里是特征长度，界面面积/某一相占的体积，单位是m^-1
        const dimensionedScalar Nl
        (
            gSum((areaDensity*mesh.V())())
           /(
               gSum
               (
                   ((areaDensity*from)*mesh.V())()
               )
             + dimensionedScalar("SMALL", dimless, VSMALL)
            )
        );


        if (mesh.time().outputTime() && debug)
        {
            areaDensity.write();
            Tmask.write();
            volScalarField mKGasDot
            (
                "mKGasDot",
                massFluxEvap*areaDensity*Nl*from
            );
            mKGasDot.write();
        }
			//返回可以用在计算模型中的相变速率，kg/m^3 s
        return massFluxEvap*areaDensity*Nl*from;
    }
    else
    {
        return tmp<volScalarField> ();
    }
}
```

总结下来这里有几处处理值得借鉴：

1. 使用mask将相变地点发生在界面上。
2. 用了Nl来使模型输出应用在模型中。（常见做法）
3. 用了大量的判断来保证相变的方向正确。

4. 将相方程RHS处理成源项的格式

下面分项分析：

1. `forAll(List, i)`第一个参数是一个`List`，温度`volScalarField也`是一个`List`，包含了所有的体心网格列表，`i`表示每次循环中读取到的网格。这里是对各相的导数进行判断，乘积不为0即说明在界面区域，然后判断两相分数之和是否到达了设定的阈值，随后将`Tmask`设为1，之后在讲这个mask与`T-Tsat`相乘，达到了限定相变位置的目的。但是`iCEF`中使用的`interfaceProperties`提供了`nearInterface()`的方法，虽然在相变模型中不好直接访问，但可以直接用alpha场构建类似的mask。

   > 经过测试，这里的`Tmask`将相变区域仅限值在气液相变区域，这导致纯液或纯气区域一直不会发生相变，例如加热平板上一直不会有气泡，这里考虑将这个限制去掉

```c++
	Foam::tmp<Foam::volScalarField>
	Foam::interfaceProperties::nearInterface() const
	{
	return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}
```

2. 这个式子从原理上是正确的，见iCEF分析，只是这个表达式可以有更多的改进。

 详见论文“An algorithm to calculate interfacial area for multiphase mass transfer through the volume-of-fluid method ”

3. 因为这里涉及到多相之间的相变，而本例不需要考虑这么多，应用时根据情况判断是否要考虑的这么仔细。
4. `iCEF`中处理成`mdotAlphal`的方式确实很蠢，如果代码变动方便的话应该考虑处理成`Su`，`Sp`的格式。这里对于压力方程中不需做太多更改，仅需改变相方程，然后将`Alhpal`结尾的方法改成使用源项。



