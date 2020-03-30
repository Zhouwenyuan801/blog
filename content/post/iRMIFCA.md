---

title: "IRMIF中接触角的设置"
date: 2020-03-30T13:03:32+08:00
lastmod: 2020-03-30T13:03:32+08:00
draft: false
keywords: ["Contact Angle","OpenFOAM"]
description: "icoReactingMultiphaseInterFoam中接触角的设置"
tags: ["Contact Angle","OpenFOAM","code"]
categories: ["OpenFOAM"]
author: "ZZQ"
---

<!--more-->

## 前言

在使用`icoReactingMultiphaseInterFoam`计算微重力下边界浸润情况时，原例中使用的是`zeroGradient`，并没有形成接触角的条件，因此改而采用`constantAlphaConotactAngle`边界条件。但是在应用过程中出现了明显的问题，即无论如何修改各个组分的边界条件，或者更改组分的类型，虽然能够顺利进行就算，但是页面边界始终未出现应有的角度。而参考`interFoam`给出的例程`capilaryFlow`，给定接触角后`cACA`边界是可以进行接触角修正的。因此，这里详细看一下这种边界处理接触角的过程，以及为什么求解器中不能正确还原接触角。



关于OpenFOAM中边界条件的处理可以参考这两篇博客[[1]](http://xiaopingqiu.github.io/2016/04/02/Boundary-conditions-in-OpenFOAM1/)[[2]](https://www.cfd-online.com/Forums/openfoam-programming-development/129271-how-boundary-conditions-called-openfoam-solvers.html)

大体上，边界条件是在物理场读入初始化时同时被创建并修正，在每一次涉及到物理场方程的运算时，也被自动修正。

>Those BC types that need updating and evaluation at each time step, they either have updateCoeffs() or evaluate() in the class definition, virtual of course.
>
>When a fvMatrix<type> Eqn is created, the updateCoeffs() is called in its constructor. Then in Eqn().solve() the function correctBoundaryConditions() is called, which contains evaluate().
>
>Now the two functions are implemented differently in every BC, and the solver would also execute different things out according to each BC type defined.

边界条件需确定两件事：1.物理场在边界上的取值。2.物理场在边界的变化情况（梯度）。这两者分别在对流项和扩散项的离散过程中被用到：
$$
\int_v\nabla \cdot (\rho U\phi)dV = \sum_f m_f \phi_f
$$

$$
\int_v \nabla \cdot(\Gamma \nabla \phi)dV  = \sum_f (\Gamma \nabla \phi_f)\cdot S_f 
$$

其中，$\phi_f$是面上的值，$S_f$是面向量，为面法向方向，大小为面积的矢量。当网格面包含物理边界时，就需要用到边界条件。

在OpenFOAM中，边界上的值、梯度通过统一的表达式表达：
$$
\phi_f = A_1 \phi_c + B_1
$$

$$
\nabla \phi_f = A_2 \phi_c + B_2
$$

其中，$\phi_c$是边界相邻网格中心的值，$A_1,A_2,B_1,B_2$是边界类中定义的系数。在类中名字分别为`valueInternalCoeffs `,`valueBoundaryCoeffs `,`gradientInternalCoeffs `,`gradientBoundaryCoeffs `。不同的边界条件，四个系数的表达式也不同。

通过调用`updateCoeffs`和`evaluate`函数完成边界值的更新，作者给出的解释

>\- on correctBoundaryConditions() for a field
>
>\- on updateCoeffs() at matrix creation
>
>correctBoundaryConditions is also called after the linear solver call automatically.

`correctBoundaryConditions()`中调用了了`evaluate`。一般来说，复杂边界多使用`updateCoeffs`，简单边界用`evaluate`。

OpenFOAM中常见的边界条件多是由几类基本条件衍生的，最常见的是三类边界条件：`Dirichlet`，`Neumann`和混合边界条件。

对于`Dirichlet`条件，即边界上给定数值，
$$
x_p = a
$$

$$
\nabla x_p = -\Delta \cdot x_c + \Delta \cdot a = (C-x_c)\cdot \Delta
$$

其中，$x_p$是边界值，$x_c$是相邻网格体心值，$\Delta$是网格中心至边界的倒数，$a$为常数。

对应的，四个常数的取值为
$$
A_1 = 0; A_2 = C; B_1 = -\Delta; B_2 = \Delta\cdot a
$$
在`fixedValue`边界的代码中，

```c++
template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}

template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}

template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}
```

这里出现的`*this`指边界本身，即边界上物理量的取值。

而`zeroGradient`边界属于`Neumann`边界，即给定边界梯度为0

$$
x_p = x_c
$$

$$
\nabla x_p = 0
$$

对应的，
$$
A_1 = 1; A_2 = 0; B_1 = 0; B_2 = 0
$$
在实现时，代码为，

```c++
template<class Type>
tmp<Field<Type> > zeroGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );
}

template<class Type>
tmp<Field<Type> > zeroGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

template<class Type>
tmp<Field<Type> > zeroGradientFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

template<class Type>
tmp<Field<Type> > zeroGradientFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}
```

另一类`Neumann`边界为你`fixedGradient`，即给定边界值，

$$
x_p = x_c + \frac{G}{\Delta}
$$

$$
\nabla x_p = G
$$

即给定边界处梯度。

另一类混合边界为上述两种的结合。除此之外，大部分边界都是基于这三种发展的。

本例中拟使用的`constantAlphaConotactAngle`边界就是基于`fixedGradient`边界演化的。下面具体看改边界的头文件。

```c++
class constantAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;
 		...        
        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;
}
Foam::tmp<Foam::scalarField>
Foam::constantAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return tmp<scalarField>::New(size(), theta0_);
}
```



除了构造函数外，该类的唯一作用就是定义接触角`theta`的计算方法，这里是读取给定的数值并返回常量。

其父类`alphaContactAngleTwoPhaseFvPatchScalarField`是更一般的方法，其中定义了大部分功能，只是把接触角计算留为纯虚，除常数子类外，还有动态接触角、温度相关接触角等其他派生类，可在运行时根据给定关键字动态选取。

父类的构造为

```c++
class alphaContactAngleTwoPhaseFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
public:

    // Abstract class, no runtime information

    //- Alpha limit options
    enum limitControls
    {
        lcNone,
        lcGradient,
        lcZeroGradient,
        lcAlpha
    };
    static const Enum<limitControls> limitControlNames_;
    limitControls limit_;
...
    // Member functions

        //- Return the contact angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const = 0;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );
}
```

使用该类，要规定边界的限制方法，常用的就是使用`Gradient`方法限制。另外，定义了`evaluate`方法，前面提到过，这就是修正边界场时使用的函数。

```c++
void Foam::alphaContactAngleTwoPhaseFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (limit_ == lcGradient)
    {
        gradient() =
        patch().deltaCoeffs()
       *(
           max(min
           (
               *this + gradient()/patch().deltaCoeffs(),
               scalar(1)), scalar(0)
           ) - *this
       );
    }
    else if (limit_ == lcZeroGradient)
    {
        gradient() = 0.0;
    }

    fixedGradientFvPatchScalarField::evaluate();

    if (limit_ == lcAlpha)
    {
        scalarField::operator=(max(min(*this, scalar(1)), scalar(0)));
    }
}
```

如代码所述，若选定`Gradient`方法，该函数起到的作用是对边界的梯度值更新：
$$
\nabla \alpha' = \Delta ([\alpha_p +\nabla \alpha/\Delta ]_{0-1} -\alpha_p)
$$
即将边界相分数限制在合理范围内。阅读求解器代码，这个函数在构造场和`alpha1.Eqn`求解结束后被调用。然而这并不是接触角的修正过程，更没有出现`theta()`和界面形状的关联。因此，需要重新寻找调用该边界`theta()`函数的地方。

经查，对接触角的修正出现在`interfaceProperties.C`中`correctContactAngle`函数中。要理解代码，先要明白是如何通过接触角修正界面形状的，这里可以参见[Grunding的文章](http://arxiv.org/abs/1907.05054)和[这篇博客](https://www.cfd-online.com/Forums/openfoam-solving/81101-interfoam-contact-angle.html)。

在使用VOF的求解器中，接触角在计算界面曲率时引入，具体的方法是更正壁面法向和界面法向的夹角
$$
N_y / \sqrt{N_x^2 + N_z^2} = \cos(\theta)
$$
或
$$
\vec{n_{f}}\cdot\vec{ n_w }= \cos(\theta)
$$
其中，$\theta$是给定的接触角，$\vec{n_{f}}$是界面法向，$\vec{ n_w }$是壁面法向。那么求解器中应该是先计算出流场内的界面法向，再根据边界条件对边界值调用`theta()`修正，这也是`icoRMF`中最可能出错的地方。那么下面具体看`interFoam`是如何修正的。


```c++
void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();
```
传递进函数的是界面法向的边界场，以及相分数梯度的边界场。随后函数调用了相分数的边界场，以及目前网格的所有边界。

```c++
    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );
            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle
```
对于相分数的每一个边界，判断它是不是`alphaContactAngleTwoPhaseFvPatchScalarField`或其子类，如果是的话，获得明确指向该边界的引用，关于`refCast`其实是`dynamic_cast`，和`const_cast`的解释见[[3]](https://en.cppreference.com/w/cpp/language/dynamic_cast)[[4]](https://stackoverflow.com/questions/13783312/how-does-dynamic-cast-work)[[5]](https://en.cppreference.com/w/cpp/language/const_cast)[[6]](https://stackoverflow.com/questions/19554841/how-to-use-const-cast/19554871)。随后拿到接触角，注意，这里传递的是速度的边界场（后面试过其他物理场都不行）和界面法向在某个边界上的值，之后转化为弧度。随后再拿到边界上壁面的法向矢量。
```c++
            const scalarField a12(nHatp & nf); //其实是cos beta，beta是当前的接触角
            const scalarField b1(cos(theta)); //cos theta 给定的接触角

            scalarField b2(nHatp.size());  
            forAll(b2, facei)  // cos (beta - theta)，应修改的角度
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}
```
随后构造的变量如注释表示，这里其实是几何问题，即已知当前、预定接触角，以及目前界面法向和边界法向矢量，求修正后的界面矢量，画图表示为，

![image-20200330140616397](https://tva1.sinaimg.cn/large/00831rSTly1gdbx6nqxkej30og0cmgnr.jpg)

然后将梯度更新再确保有界。

更新接触角是在计算曲率的`calculateK`函数中调用的，该函数中，计算界面法向场后即调用接触角函数，之后用更新后的法向场计算曲率。而`calculateK`在构造`interfaceProperties`对象时调用，在每次`correct`时也调用你。在`interFoam`中，所创建的`mixture`继承自`interfaceProperties`，因此，在`createFields.H`中，求解相方程后均修正了界面形状。

这套处理逻辑应该是有效的，如前所述，`interFoam`可以计算毛细流。那么就要看一下`iRMIF`中是哪些环节出了问题。

在本例中，求解器并没有创建`mixture`，而是创建了`multiphaseSystem`对象`fluid`，对象中不包含`interfaceProperties`中的方法，而本身的`correct`方法（在`T.Eqn.`中使用）作用时更新温度场。但是在你`phaseSystem`中重新写了`nHatfv`函数以计算界面法向。

```c++
Foam::tmp<Foam::surfaceVectorField> Foam::phaseSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}
```

可见，对比`interFoam`，明显缺少了边界修正的过程。改进方法有两种：1. 构造`interfaceProperties`对象，使用标准接口处理接触角。2.修改`nHatfv`函数，植入接触角修正过程。初步评估第二种对体系改动小一些，因此选择第二种，增加如下代码，即可实现接触角的修正

```c++
    tmp<surfaceVectorField> tnhatfv = gradAlphaf/(mag(gradAlphaf) + deltaN);
    // Face unit interface normal
    surfaceVectorField& nhatfv = tnhatfv.ref();
    surfaceVectorField::Boundary& nHatb = nhatfv.boundaryFieldRef();
    const surfaceVectorField::Boundary& gradAlphafb = gradAlphaf.boundaryField();

    const volScalarField::Boundary& a1bf = alpha1.boundaryField();
    // const volScalarField::Boundary& a2bf = alpha2.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();


    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(a1bf[patchi]))
        {

            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        a1bf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];

            volVectorField U_(
                IOobject
                (
                    "U_",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
            );
            forAllConstIters(phaseModels_, iter)
            {
                U_ += iter()() * iter()->U();
            }

            const scalarField theta
            (
                degToRad(acap.theta(U_.boundaryField()[patchi], nHatp))
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }
            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphafb[patchi]);
            acap.evaluate();
        }
```

值得注意的一点是，该函数返回的是`tmp<surfaceVectorField>`。这就导致一个问题，即使在原函数不变的基础上，定义一个中间变量`surfaceVectorField`，再返回这个场，编译可以通过，但运行时会报错。因此，只能先构建一个`tmp`封装的对象，再对它的引用操作。关于`tmp`机制可以参见[这篇博客](https://marinecfd.xyz/post/openfoam-tmp/)。

