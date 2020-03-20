---
title: "IRMIF中的多流体模型"
date: 2020-03-04T15:34:30+08:00
lastmod: 2020-03-04T15:34:30+08:00
draft: false
keywords: ["E&C","OpenFOAM"]
description: "IcoReactingMultiphaseInterFoam中多流体模型结构"
tags: ["E&C","OpenFOAM","code"]
categories: ["OpenFOAM"]
author: "ZZQ"
---

<!--more-->

`interCondensatingEvaporatingFoam`在计算低温流体时温度场发散，在界面区域温度出现小范围波动后，迅速在气液界面区域出现极高温和极低温，高至几百K，低至负几千K（最早时温度场在迭代几步之后出现发散，之后偶然有一次可以计算，最近再试依旧发散，所以最终归结为其对于低温流体温度、物性的计算不稳定）。初步归纳发散的原因有以下几个可能，

1. 温度控制方程太简陋，仅仅是Laplace方程；
2. 质量方程中源项处理不准确，未考虑源项为负时显式处理。

那么最终还是要回到使用`icoReactingMultiphaseInterFoam`上，这个求解器有以下几个特点：

1. 对热物性方面处理更佳缜密，沿用了`OpenFOAM`的热物性库；
2. 对多相（大于等于2）的处理更精确，包括相内多组分的求解；
3. 使用液氧计算简单沸腾现象通过；
4. 模型结构更复杂，内容更多。

之前已经分析过该求解器中的相变模型，但要更准确的使用，还需要再探讨以下这些内容：

1. 求解器计算过程
2. 多流体模型
3. 相模型

对于辐射模型和化学反应过程这里不关心。这篇争取解释好前两点。

## 求解器脉络

先按初始化各个物理场，`createFields.H`中创建的场包括：压力、速度、温度，`rhoCp`，`K`，密度，`rhoPhi`。创建的指针包括：多流体模型，湍流模型，辐射模型的指针。

关于时间步长的确定有两个准则，一是常规采用的有限体积法下计算的库朗数[[1]](https://marinecfd.xyz/post/courant-number-in-unstructured-grids/)

```c++
    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
```

还引用了`alphaCourantNo.H`计算相界面处的库朗数，

```c++
    scalarField sumPhi
    (
        fluid.nearInterface()().internalField()
       *fvc::surfaceSum(mag(phi))().internalField()
    );

    alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanAlphaCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    ddtAlphaNum = fluid.ddtAlphaMax().value()*runTime.deltaTValue();

    DiNum = fluid.maxDiffNo();
```

`nearInterface()`方法返回的场在界面区域为1，其他地方为0，这里起到mask的作用。同时获得多相体系中的最大相变量`ddtAlphaNum`和扩散系数`maxDiffNo`，体系的方法分别调用体系内所存流体模型的方法，并输出最大值。

做完准备工作后，求解器分别求解`UEqn`,`PEqn.H`,`TEqn.H`和`YEqn.H`。前三个比较常见，这里主要分析耦合相变模型后方程的变化，最后的组分方程在求解相内组分迁移时引用。下面按照方程求解顺序介绍。

### `fluid.solve()`

在求解控制方程前先调用了多流体系统的`fluid.solve()`方法， 该方法定义于`multiphaseSystem`类中。该方法对于每一相采用类似`interfoam`中求解两相分布的方法，定义相迁移方程，并使用MULES修正。

同时，在求解过程中，调用`calculateSuSp()`方法，更新相变源项。源项包括`Su`，`Sp`项，其中`Su`显式定义在相方程中，`Sp`是带有$\alpha$ 的源项，通过调用`fvm::Sp`讲源项更新在系数矩阵对角线上求解。在本例中程序会判断源项的正负，若`Sp`为正值，则全部更新到`Su`中，若为负，则隐式求解，这样可以保证系数矩阵的对角占优，增强数值稳定性，更一般的方法还有调用`fvm::SuSP`，可以自动处理源项[[2]](https://marinecfd.xyz/post/openfoam-source-term-treatment/)。

### UEqn.H

求解步骤与一般的UEqn并无明显区别，关于表面张力模型后面再看。在建立UEqn并松弛之后，多了一项跟多孔介质有关的描述`fluid.addInterfacePorosity(UEqn);`。该方法在UEqn矩阵对角线上增加多孔介质模型的`S()`方法，因为不懂，所以这里不展开描述。

### YEqn.H

组分方程是求解器中重要的内容。方程描述的是单相内的组分迁移规律，各组分之间是互融的，组分迁移是自由扩散迁移。

进到文件里面先套了一层循环，冒号指表示遍历可iterate的变量。

```c++
for (phaseModel& phase : fluid.phases()):

const Foam::UPtrList<Foam::phaseModel>& Foam::multiphaseSystem::phases() const
{
return phases_;
}
```

其中，`phases_`是指针列表，指向多流体中的各个流体模型。进入循环后，在单相内构建组分源项，并调用`fluid.massSpeciesTransfer()`方法计算相间的质量变化。最后调用相模型的`phase.solveYi()`方法求解组分迁移。

```c++
forAll(Y, i)
{
  // Calculate mass exchange for species consistent with
  // alpha's source terms.
  fluid.massSpeciesTransfer(phase, Sus[i], Sps[i], Y[i].name());
}
phase.solveYi(Sus, Sps);
```

`massSpeciesTransfer()`在`MassTransferPhaseSystem`中定义，代码及注释如下

```c++
template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massSpeciesTransfer
(
    const phaseModel& phase,
    volScalarField::Internal& Su,
    volScalarField::Internal& Sp,
    const word speciesName
)
{
    // Fill the volumetric mass transfer for species
    forAllConstIters(massTransferModels_, iter)
    {
        if (iter()->transferSpecie() == speciesName)
        {
            // Explicit source，将相的变化率读取至组分，并显式求解。
            Su +=
                  this->Su()[phase.name()]
                + this->Sp()[phase.name()]*phase.oldTime();
        }
    }
}
```

`solveYi()`方法在各自的相模型中定义，对于纯质，如`purePhaseModel`，这个函数什么都不做，对于混合物，如`multiComponentPhaseModel`，这个函数在对组分场进行一些处理后，单独求解各个组分的方程，并做了MULES修正。

### TEqn.H

温度方程的形式大体是焓方程，但是和`OpenFOAM`官方给出的形式稍有差别，数学表达式可参考[这篇帖子](https://www.cfd-online.com/Forums/openfoam-solving/150535-evapvofhardt-discussion-come-join.html)，以及[Hardt](https://linkinghub.elsevier.com/retrieve/pii/S0021999108001228)的文章。温度方程中包含`fluid.heatTransfer(T)`，表证的是因传质引入的源项，该源项对每一对相变物理场更新源项为
$$
S_m = [h_{lg}-(c_{p,l}-c_{p,g})T_{std}]\dot{m}
$$


其中，$T_{std}$是`constant::standard::Tstd`，表示单位温度。考虑了这部分相变质量变化时潜热和显热的部分。

在求解温度方程后，重新调用了`fluid.correct()`方法调用每一个相的修正函数。

### PEqn.H

与`iCEF`中类似，在每一对发生相变的相方程右侧增加了传质引入的源项，

```c++
p_rghEqn +=
dmdtNet*
(
- fluid.coeffs(phase1.name())
+ fluid.coeffs(phase2.name())
);
```

`coeffs`返回的是密度的倒数，对应的方程为，
$$
\sum_{f\in\partial \Omega_i}\left ( \frac{1}{A_P}\right )\left ( \nabla^\perp_f p^{m+1}_d\right ) |S_f| - \sum_{f\in\partial \Omega_i} \phi^r_f = \frac{\dot{m}}{\rho_l}- \frac{\dot{m}}{\rho_v}
$$
在PIMPLE循环结束后，更新体系内的密度，然后将计算结果输出。

## 多流体模型

求解器中多流体模型主要由`phaseModel`，`phaseSystem`，`multiphaseSystem`，`MassTransferPhaseSystem`构成，其逻辑继承关系基本也是从前到后依次扩大。`phaseModel`是相模型，其本身继承自`volScalarField`，是表示相分数的场。同时，类中定义了大量与相绑定的物性参数计算方法，如密度、粘度等。具体有关相模型的内容后面再介绍，本节主要介绍求解器的多流体模型是如何搭建的。

在初始化物理场时，创建的是指向`multiphaseSystem`的指针，

```c++
Info<< "Creating multiphaseSystem\n" << endl;
autoPtr<multiphaseSystem> fluidPtr = multiphaseSystem::New(mesh);
multiphaseSystem& fluid = fluidPtr();
```

该类继承自`phaseSystem`，这两个类的区别如下，



|          |                 multiPhaseSystem                 |                         phaseSystem                          |
| :------: | :----------------------------------------------: | :----------------------------------------------------------: |
|   父类   |                   phaseSystem                    |         basicThermo, <br>compressibleTransportModel          |
| 包含的表 | SuSpTable<br>scalarTable<br>compressionFluxTable | phasePairTable<br>phaseModelTable<br>sufaceTensionTable<br>interfacePorousTable |

从功能上讲，`phaseSystem`基本已经构建完成了多流体系统的框架，包括了体系内相对，相模型和子模型表。而`multiphaseSystem`虽然在名字上多了`multi`，但是增加的内容大多是与传热、求解相方程有关的，包括前述solve()`方法。

在`createField.H`中，使用的是`New`(selector)完成创建指针，那么得到的是RTS机制获得的指针。以`evaporatingMultiComponent`案例为例，在`constant/phaseProperties`中`type`关键字为`massTransferMultiphaseSystem`，这说明后者是前者的派生类。该类中包含了相变模型(`interfaceCompositionModel`)及具体计算传热传质量的表达式。

当构造一个`MassTransferPhaseSystem`类时，构造的顺序是`phaseSystem`->`multiphaseSystem`->`MassTransferPhaseSystem`。下面按照构造顺序依次分析多流体系统的构建过程。

### phaseSystem的构造



构造时先调用其父类`basicThermo`的构造函数，保存`mesh`的引用，并在`constant/phaseProperties`中查找`phases`字段并存储至`phaseNames`列表中，同时生成压力、温度场、相分数场。

完成`basicThermo`构建后，调用`generatePhaseModels`初始化`phaseModels_`这张hash表。这个方法需要的参数是`wordList& phaseNames`

```c++
    phaseModelTable phaseModels;

    for (const word& phaseName : phaseNames) //对读到的每一相都创建模型
    {
        phaseModels.insert
        (
            phaseName,
            phaseModel::New //phaseModel中New时需要声明所属的phaseSystem
            (
                *this,
                phaseName
            )
        );
    }

    return phaseModels;
```

`phaseModel`的创建在后面讨论。这样在`phaseModels_`中初始化了一张`key`为相名称，`val`为指针的表。

之后，寻找`phaseProperties`中是否有`surfaceTension`或者`interfacePorous`子字典，如果有的话，调用`generatePairsAndSubModels`创建相对和子模型。

```c++
    if (found("surfaceTension"))
    {
        generatePairsAndSubModels
        (
            "surfaceTension",
            mesh_,
            surfaceTensionModels_
        );
    }
```

其中，表面张力模型表的类型为

```c++
        typedef
            HashTable
            <
                autoPtr<surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;
```

对应的`generatePairsAndSubModels`方法定义为

```c++
    dictTable modelDicts(lookup(modelName));

    generatePairs(modelDicts);

    createSubModels(modelDicts, models);
```

首先在文件中搜索`modlename`，比如，在这里是`surfaceTension`，例程中为，

```c++
surfaceTension
(
    (gas and liquid)
    {
        type            constant;
        sigma           0.07;
    }
);
```

然后将该段feed给`HashTable<dictionary, phasePairKey, phasePairKey::hash> modelDicts`。该表是由`phasePairKey`和`dictionary`组成的Hash表。由例程可见，`gas and liquid`用来初始化一个`phasePairKey`，而后面的内容作为一个字典存至表中。

---

这里额外分析一下`phasePairKey`，`phasePair`，和`orderedPhasePair`。继承关系为从前到后。

`phasePairKey`包含了一对字符，和表征字符间是否有顺序的布尔值。重定义了IO和比较逻辑，例如，当读入的是`xx and xx`时，`bool = false`；当读入`xx to xx`时，`bool = ture`。同时，还定义了一个`hash()`，用以生成存储的字符对对应的hash值。`phasePair`继承自`phasePairKey`，增加了两个`phaseModel`的引用和接口，但是默认将`ordered`设置为`false`。自然地，`orderedPhasePair`描述的是两相间顺序有意义的相对。

---



之后调用`generatePairs`创建相对，

```c++
    forAllConstIters(modelDicts, iter)//按Key生成Pair，并储存
    {
        const phasePairKey& key = iter.key();//获得当前循环的Key

        // pair already exists
        if (phasePairs_.found(key))
        {
            // do nothing ...
        }
        else if (key.ordered())
        {
            // New ordered pair
            phasePairs_.insert //该表是phasePairKey和phasePair组成的
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
        else
        {
            // New unordered pair
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
```

生成`phasePair`后，填补表征子模型的表（比如哪个相对应采用哪个表面张力模型），调用的是`createSubModels`，

```c++
template<class modelType>
void Foam::phaseSystem::createSubModels
(
    const dictTable& modelDicts,
    HashTable
    <
        autoPtr<modelType>, //表面张力: interfaceCompositionModel
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    forAllConstIters(modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        models.insert
        (
            key,
            modelType::New
            (
                iter.val(),  //按照读进来的字典文件创建表面张力模型
                phasePairs_[key] //key对应的相对
            )
        );
    }
}
```



关于表面张力模型的创建后面再详细介绍。

子模型创建完成后，调用`generatePairsTable`完成体系内相对的构建，以下方法保证创建出的`totalPhasePair`表中元素互异。

```c++
forAllConstIters(phaseModels_, phaseIter1)//从相模型字典中读取相1
    {
        forAllConstIters(phaseModels_, phaseIter2)//读取相2
        {
            if (phaseIter1()->name() != phaseIter2()->name())//判断不同相
            {
                phasePairKey key //构建有1至2的有序Key
                (
                    phaseIter1()->name(),
                    phaseIter2()->name(),
                    true
                );

                phasePairKey keyInverse //构建反向Key
                (
                    phaseIter2()->name(),
                    phaseIter1()->name(),
                    true
                );

                if
                (
                    !totalPhasePairs_.found(key)
                 && !totalPhasePairs_.found(keyInverse)
                ) //若正向反向key均为包含在total中，构建phasepair
                {
                    totalPhasePairs_.set
                    (
                        key,
                        autoPtr<phasePair>
                        (
                            new phasePair
                            (
                                phaseModels_[key.first()],
                                phaseModels_[key.second()]
                            )
                        )
                    );
                }
            }
        }
    }
```

需要注意的是，`phasePairs_`是在`generatePairs`方法中使用，是靠读取字典文件中的关键字一一创建。而`totalPhasePairs`则是按照排列组合把体系包含的所有相对都创建一遍。

至此，就完成了`phaseSystem`的构建。该类可提供多流体系统中物性的接口，比如密度，

```c++
Foam::tmp<Foam::scalarField> Foam::phaseSystem::rho(const label patchI) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpRho
    (
        iter()().boundaryField()[patchI]
      * iter()->rho()().boundaryField()[patchI]
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpRho.ref() +=   //加权求和体系内各model的密度
        (
            iter()().boundaryField()[patchI]
          * iter()->rho()().boundaryField()[patchI]
        );
    }

    return tmpRho;
}
```

同时，有关表面张力的计算也在`phaseSystem`中完成，比如计算曲率的`K()`，以及计算表面张力的`surfaceTension()`。定义的表面张力模型只是返回一个表面张力系数`sigma_`

### multiphaseSystem的构建

`multiphaseSystem`在调用`phaseSystem`构造函数的基础上，读取了一些MULES的控制参数，同时初始化了`Su`，`Sp`两张源项表。

该类额外提供的功能包括

- 求解体系内相迁移方程
- 计算相方程中源项的方法
- 以`phase[i]`方法获取`phaseModel`或者体积分数场

### MassTransferPhaseModel

该类是继承自`multiphaseSystem`的RTS类，通过`multiphaseSystems.C`加入RTS表

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

在构造时，先调用父类的构造函数，同时，利用前面提到的创建子模型的方法`generatePairsAndSubModels`创建传质模型表，

```c++
this->generatePairsAndSubModels("massTransferModel", massTransferModels_);
```

以例程为例，字典如下所示

```c++
massTransferModel
(
    (liquid to gas)
    {
        type            Lee;
        species         vapour.gas;
        C               8;
        alphaMin        0.0;
        alphaMax        0.2;
        Tactivate       90.147;
    }
);
```

构造完成后，声明`dmdt`（传质速率）场，并且声明需要保存。

```c++
    forAllConstIters(massTransferModels_, iterModel)
    {
        if (!dmdt_.found(iterModel()->pair()))
        {
            dmdt_.set
            (
                iterModel()->pair(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("dmdt",iterModel()->pair().name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                )
            );
        }
    }
```

该类中包含了全部的与计算传热相关的函数，比如计算潜热的`calculateL`，计算传热并切更新源项的`heatTransfer`，以及计算组分间传质的`massSpeciesTransfer`。



至此，多流体模型的框架久搭建完成了。

## PhaseModel构建

`phaseModel`这个类在多流体模型中起到中坚的作用。多流体模型只是将各个相模型的物性加权求和后返回，真正确定流体物性的是`phaseModel`类。

前文中提到，`phaseModel`是在`phaseSystem`构造函数中调用`generatePhaseModel`完成的，调用的是New函数。

```c++
Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::New
(
    const phaseSystem& fluid,
    const word& phaseName
)
{
    const dictionary& dict = fluid.subDict(phaseName); //读取相名对应的子字典

    const word modelType(dict.get<word>("type")); //获得给定的相模型名称

    Info<< "Selecting phaseModel for "
        << phaseName << ": " << modelType << endl;

    const auto cstrIter = phaseSystemConstructorTablePtr_->cfind(modelType); //在RTS表中寻找对应类的New指针

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "phaseModel",
            modelType,
            *phaseSystemConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return cstrIter()(fluid, phaseName);//返回相模型指针
}
```

例程中给出的形式为

```c++
liquid
{
    type            pureMovingPhaseModel;
}

gas
{
    type            multiComponentMovingPhaseModel;
}
```

这两个类别是在`phaseModels.C`中定义的，

```c++
makePhaseTypes
(
    MovingPhaseModel,
    PurePhaseModel,
    phaseModel,
    rhoThermo,
    pureMovingPhaseModel // Name of the phase type
);
makePhaseTypes
(
    MovingPhaseModel,
    MultiComponentPhaseModel,
    phaseModel,
    rhoReactionThermo,
    multiComponentMovingPhaseModel
);
```

`makePhaseTypes`是在`makePhaseTypes.H`中定义的宏，

```c++
#define makePhaseTypes(MomemtumType, CompType, Phase, Thermo, Name)         \
                                                                            \
    namespace Foam                                                          \
    {                                                                       \
        typedef Foam::MomemtumType                                          \
        <                                                                   \
            Foam::CompType                                                  \
            <                                                               \
                Foam::Phase,                                                \
                Foam::Thermo                                                \
            >                                                               \
        >                                                                   \
        Name;                                                               \
                                                                            \
        addNamedToRunTimeSelectionTable                                     \
        (                                                                   \
            phaseModel,                                                     \
            Name,                                                           \
            phaseSystem,                                                    \
            Name                                                            \
        );                                                                  \
    }
```

该宏声明特定类的继承关系，并加入RTS表中。以`multiComponentMovingPhaseModel`为例，其完整的类型应该是`MovingPhaseModel<MultiComponentPhaseModel<phaseModel, rhoReactionThermo>> multiComponentMovingPhaseModel` 。

从内向外分析，`phaseModel`是上述模型基类，通过读取时间步文件夹内的`alpha.(phaseName)`初始化。该类中定义了几乎所有的物性参数，包括`alpha`（热扩散系数），`rho`（密度），`Cp`（定压比热），`Cv`（定容比热），`Cpv`，`CpByCpv`，`nu`（粘度），`Y`（组分场），`phi`（体积通量），`alphaPhi`（相体积通量），`U`。但是全部的物性都是通过`thermo()`获得的，比如比热，`return thermo().Cp();`。而`thermo()`是纯虚函数，

```c++
virtual const rhoThermo& thermo() const = 0;
```

因此，`phaseModel`类虽定义了相分数场，但是物性需具体实例化后才能获取。

`rhoReactionThermo`继承自`rhoThermo`，定义了纯虚函数`composition()`，用以返回混合物中的组分。在构造时，调用的是`basicThermo::New<rhoReactionThermo>(mesh, phaseName, dictName)`。

`rhoThermo`继承自`fluidThermo`，定义了密度、粘度和压缩性的成员变量，并给出了数据接口。

`fluidThermo`继承自`basicThermo`和`compressibleTranspooortModel`，给出了粘度的定义。

`basicThermo`继承自`IOdictionary`，存储了压力，温度的引用，同时定义了热扩散系数，相名称，以及包括比热等纯虚的数据接口。除此之外，还定义了一些返回字符的函数，例如

```c++
static word phasePropertyName
  (
  const word& name,
  const word& phaseName
)
{
  return IOobject::groupName(name, phaseName);
}

word phasePropertyName(const word& name) const
{
  return basicThermo::phasePropertyName(name, phaseName_);
}
const Foam::word Foam::basicThermo::dictName("thermophysicalProperties");
```

这些是为了后续方便。

那么回头看一下构造时调用的New函数，该方法定义在`basicThermoTemplate.C`中，

```c++
        phaseThermo::New //MultiComponentPhaseModel中调用的构造函数
        (
            fluid.mesh(),
            phaseName,
            basicThermo::phasePropertyName(basicThermo::dictName, phaseName)
        ).ptr()
 
template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
  	IOdictionary thermoDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    auto cstrIter =
        lookupThermo<Thermo, typename Thermo::fvMeshDictPhaseConstructorTable>
        (
            thermoDict,
            Thermo::fvMeshDictPhaseConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, phaseName, dictName));
}
```

其中参数`basicThermo::phasePropertyName(basicThermo::dictName, phaseName)`代入后是`thermophysicalProperties.gas/liquid`位于`constant`文件夹下。使用`lookupThermo`方法获得指针，

```c++
     if (thermoDict.isDict("thermoType")) //判断文件中是否存在thermoType字段
     {
         const dictionary& thermoTypeDict = thermoDict.subDict("thermoType");
  
         Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;
  
         if (thermoTypeDict.found("properties")) //如果有properties关键字
         {
             std::initializer_list<const char*> cmptNames
             {
                 "type",
                 "mixture",
                 "properties",
                 "energy"
             };
  
             // Construct the name of the thermo package from the components
             const word thermoTypeName
             (
                 thermoTypeDict.get<word>("type") + '<'
               + thermoTypeDict.get<word>("mixture") + '<'
               + thermoTypeDict.get<word>("properties") + ','
               + thermoTypeDict.get<word>("energy") + ">>"
             );
  
             return lookupThermo<Thermo, Table>
             (
                 thermoTypeDict,
                 tablePtr,
                 cmptNames,
                 thermoTypeName
             );
         }
         else  //若没有
         {
             std::initializer_list<const char*> cmptNames
             {
                 "type",
                 "mixture",
                 "transport",
                 "thermo",
                 "equationOfState",
                 "specie",
                 "energy"
             };
  
             // Construct the name of the thermo package from the components
             const word thermoTypeName //根据字典构建名称
             (
                 thermoTypeDict.get<word>("type") + '<'
               + thermoTypeDict.get<word>("mixture") + '<'
               + thermoTypeDict.get<word>("transport") + '<'
               + thermoTypeDict.get<word>("thermo") + '<'
               + thermoTypeDict.get<word>("equationOfState") + '<'
               + thermoTypeDict.get<word>("specie") + ">>,"
               + thermoTypeDict.get<word>("energy") + ">>>"
             );
  
             return lookupThermo<Thermo, Table> //从表中返回名称对应的指针
             (
                 thermoTypeDict,
                 tablePtr,
                 cmptNames,
                 thermoTypeName
             );
         }
     }
```

可见返回的是有`thermoType`字典内关键字规定的特定类型的指针，例程中对应，

```c++
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       const;
    thermo          hConst;
    equationOfState incompressiblePerfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}
```

关于hash表的构建，以及类的继承关系、读取指针有两篇文章写的更详细些[[3]](http://xiaopingqiu.github.io/2016/06/25/thermophysics2/)[[4]](http://xiaopingqiu.github.io/2016/06/25/thermophysics3/)。以上述字典为例，最终构建的是，

```c++
heRhoThermo
  <
  	rhoThermo,
		multiComponentMixuture,
		<
 constTransport<species::thermo<hConstThermo<imcompressiiblePerfectGas<specie>>,sensibleInternalEnergy>>
      >
  >
```

所以`phaseModel`中`thermo`指向的最终是`heRhoThermo`。`heRhoThermo`继承自`heThermo`，并以第一个模版参数为基类，第二个模版参数为`MixtureType`。到这里可能有点乱，需要明确的是，`rhoReactionThermo`调用了`basicThermo`的selector，构造了上面复杂的类，此外，还另外定义了一些方法。而其他的物性计算方法，有字段中给定的关键字确定。在给定关键字后，要在文件中给出各个物性需要的参数，以例程中气体为例，

```c++
species
(
    air
    vapour
);

inertSpecie air; //惰性组分

vapour
{
    specie
    {
        nMoles      1;
        molWeight   32;
    }
    equationOfState //因为这里选的是不可压理想气体，需要给定参考压力
    {
        pRef         1e5;
    }
    thermodynamics
    {
        Hf          0;//-1.338e7; //[J/kg]
        Cp          935.5;
    }
    transport //这里选用的是constant只需给定一个数
    {
        mu          7e-06;
        Pr          0.7;
    }
}

air
{
    specie
    {
        nMoles      1;
        molWeight   28.9;
    }
    equationOfState
    {
        pRef         1e5;
    }
    thermodynamics
    {
        Hf          0;
        Cp          900;
    }
    transport
    {
        mu          1.8e-05;
        Pr          0.7;
    }
}
```

可选的物性方法为见官方文档[[5]](https://cfd.direct/openfoam/user-guide/v6-thermophysical/)，更精确的，小范围的物性计算方法是在温度区间内使用多项式拟合，对应的需要将`transport`,`thermodynamic`,`equation of state`字段修改成`polynomial`,`hPolynomial`,`icoPolynomial`。

抄录各个方法所需参数如下

```c++
    transport
    {
        muCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
        kappaCoeffs<8>  ( 2000 -0.15 0.023 0 0 0 0 0 );
    }
```

$$
\mu  = 1000 - 0.05 T + 0.003 T^2
$$

$$
\kappa = 2000 - 0.15 T + 0.023 T^2
$$



```c++
    thermodynamics
    {
        Hf              0; //Heat of formation
        Sf              0; //Standard entropy
        CpCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
```

$$
Cp = 1000 - 0.05 T + 0.003 T^2
$$

```c++
		equationOfState
		{
        rhoCoeffs<8>    ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
```

$$
\rho = 1000 - 0.05 T + 0.003 T^2
$$

多项式拟合的优点很突出：不需要很高的阶数在给定小区间内精度就很高，缺点是在大区间内精度急剧下降。

完全创建`thermo`后，就完成了宏中`CompType`的创建，求解器中只有两类：`MultiCoponentPhaseModel`和`PurePhaseModel`分别对应混合物和纯质。这两个类进一步定义了相内组分场的计算和接口。纯质直接返回组分场`Y_`，而混合物通过`solveYi()`求解组分，前面已有分析。

定义纯质或混合物后，这个类作为模版参数传入`MomentumType`中，本例中有`MovingPhaseModel`和`StaticPhaseModel`，对应的，就是流体和固体。固体中的速度，通量都为0，而流体中的返回的是注册在`mesh`中的速度、通量等。

在`thermophysicalProperties.X`中按字段修改后并不能直接使用该物性模型，在`InterfaceCompositionModels/InterfaceCompositionModels.C`中，求解器只预定义了有限几种组合形式。如若按上述表述，应在该文件中增加引用

```c++
#include "icoPolynomial.H"
#include "polynomialTransport.H"
#include "hPolynomialThermo.H"
```

并增加新的类型的声明和定义

```c++
    typedef
        polynomialTransport
        <
            species::thermo
            <
                hPolynomialThermo
                <
                    icoPolynomial<specie>
                >,
                sensibleEnthalpy
            >
        > polyRhoHThermoPhysics;
    // Polynomial fitting
        makeInterfaceContSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polyRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );
        makeInterfaceContSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polyRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            polyRhoHThermoPhysics
        );
        makeInterfaceContSpecieMixtureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polyRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );
        makeInterfaceContSpecieMixtureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polyRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            polyRhoHThermoPhysics
        );
```

由此可见，每新增一个模型，都要在该文件中声明一遍。

至此，完整的`phaseModel`就构建完成了。