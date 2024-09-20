
# Model description

## The base model

### Dynamic Energy Budget component

`DEBBase.jl` provides a base model as ODE system. This model is a variant of DEBkiss (**citation needed**). <br>

For mathematical notation, we use a simplified DEB notation, 
where a leading capital character denotes a state variable, $k_x$ denotes a first-order 
rate canstent with respect to process $x$ and $\eta_{ab}$ denotes the efficiency 
of converting $a$ to $b$. Superscripts indicate life stage-specificity of a value 
($emb$ = embryo, $juv$=juvenile, $ad$=adult). $x_0$ indicates the baseline value of $x$ which may be affected by a stressor or other environmental or internal factors (e.g. temperature, energetic state). 
<br>
We use dot notation for derivatives over time and do not use dot notation on parameters. <br>
In line with standard DEB notation, ${x}$ denotes that $x$ is relative to a surface area and 
$[x]$ is relative to a volume (or a linearly proportional quantity, respectively). <br>
The variable naming in the source code also follows these principles as far as possible.
<br>

The ODE system has the major state variables $X$ (external resource abundance), 
$X^{emb}$ (mass of the vitellus/yolk), $S$ (structural mass), $H$ (maturity level) and $R$ (reproduction buffer). 
The intermediate variables $I$ and $A$ relate to the ammount of ingested and assimialted food, 
respectively. Maintenance costs are split into somatic maintenance $M$ and maturity maintenance $J$.
The energy flux starts with ingestion $\dot{I}$ of a resource $X$.

The mathematical model formulation follows from the basic assumptions of the standard DEBkiss model (**citation needed**):

$$

\dot{I} = \begin{cases}
    f(X)\ \{ \dot{I} \}_{max} S^{2/3}\ \ mit\ f(X) = \frac{[X]} {[X] + K_X}\ \ if\ \ X^{emb} = 0 \\
    \\
    \{ \dot{I} \}_{max}  S^{2/3}\ \ if\ \ X^{emb} > 0
\end{cases} \\
\ \\
$$

$$
\dot{X}^{emb} = \begin{cases}
    -\dot{I}\ \ if\ X^{emb} > 0 \\
    \\
    0\ \ if\ X^{emb} = 0
\end{cases}
\ \\
\ \\
\dot{X} = \begin{cases}
    -\dot{I}\ \ if\ \ X^{emb} = 0 \\
    \\
    0\ \ if\ \ X^{emb} > 0
\end{cases}
\ \\
\ \\
\dot{A} = \eta_{IA} \dot{I} \\
\ \\
\dot{M} = k_M S \\
\ \\
\dot{J} = k_J H \\
\ \\
\dot{S} = \begin{cases}
    \kappa \dot{A} - \dot{M}\ \ if\ \ \kappa \dot{A} \geq \dot{M} \\
    \ \\
    -(\dot{M} / \eta_{SA} - \kappa \dot{A})\ \ if\ \ \kappa \dot{A} < \dot{M} 

\end{cases}
\ \\
\ \\
\dot{H} = \begin{cases}
    (1 - \kappa) \dot{A} - k_J H\ \ if\ \ H < H^p \\
    \ \\
    0\ otherwise
\end{cases}
\ \\
\ \\
\dot{R} = \begin{cases}
    (1 - \kappa) \dot{A} - k_J H\ \ if\ \ H \geq H^p \\
    \ \\
    0\ \ otherwise
\end{cases}

$$

A key aspect of this model variant is that all major state variables, including structure, have the same dimension, which is assumed to represent a mass or energy pool (e.g. given in $g\ dry\ mass$, $\mu g\ C$, $J$). To retrieve a prediction of organism length, a weight-length relationship should be applied externally to the ODE. <br>
A consequence of this model formulation is the dimensionless parameter growth efficiency ($\eta_{AS}$), which replaces the volume-specific costs for growth. <br>
Surface area scaling is done through the scaling factor $S^{2/3}$, where $S$ is structural mass. 
This also leads to a different dimension of the specific ingestion rate $\{ \dot{I} \}_{max}$. 
For example, if state variables are given in $g\ dry\ mass$, 
the unit of $\{ \dot{I} \}_{max}$ is $g\ g^{2/3}\ d^{-1}$, or $\sqrt[3]{g}\ d^{-1}$. 

The zoom factor $Z$ is often used to perform corrections of DEB parameters on the basis of maximum length (**citation needed**). <br>
In line with the philosophy of using a single model currency, we define the zoom factor $Z_M$ on the basis of the model currency. 
$$
Z_M = \frac{S_{max,2}}{S_{max,1}}
$$

$S_{max}$ is maximum structural mass and can be calulated based on the equlibria of the model:

$$
S_{max} = \left( \frac{\kappa_0 \eta_{IA} \{ \dot{I} \}}{k_{M,0}} \right)^3
$$

$Z_M$ is the ratio between the maximum structural masses of two organisms, 
and is used to scale a parameter set to imply a different maximum size. 
Due to the relationship between $S_{max}$ and 
$\{\dot{I}\}_{max}$, we can establish the relationship

$$
\{ \dot{I} \}_{max,2} = \{ \dot{I} \}_{max,1}\ Z^{1/3}
$$

The zoom factor can further be applied to parameters for which it is plausible to assume that 
they scale with mass. By default, `DEBBase.jl` assumes this to be the following:


$$
\{ \dot{I} \}_{max,2}^{emb} = \{ \dot{I} \}_{max,1}^{emb}\ {Z_M}^{1/3} \\
\ \\
X^{emb}_{int,2} = X^{emb}_{int,1} Z_M \\ 
\ \\
H^p_2 = Z_M H^p_1 \\
\ \\
K_{X,2} = Z_M\ K_{X,1}
$$

`DEBBase.jl` also uses $Z_M$ to induce individual variability in DEB parameters. <br>
This is done by defining $Z_M$ as a distribution on the species level. <br>
$Z_M$ is defined as Dirac distribution if no individual 
variability should be simulated.


### Toxicokinetic-Toxicodynamic component

The base model also includes a toxicokinetic/toxicodynamic (TKTD) component. 
The TKTD component allows for an arbitray number of chemical stressors to be simulated 
as a mixture, where each stressor can act via an arbitrary combination 
of physiological modes of action (PMoA).
Toxicokinetics (TK) are implemented with a minimal model based on scaled damage (**citation needed**) 
and without interactions with body size. 
The scaled damage and associated rate parameter are specific to the PMoA, indicated by subscript $j$.
Toxicodynamics (TD) assumes a dose-response relationship between the scaled damage and the 
rate or efficiency which is affected through the corresponding PMoA. 
The intermediate variable "stress" (**citation needed**) does not explicitly appear in the model formulation. Instead, dose-response functions are directly defined to be appropriate for a given PMoA.

For PMoAs which cause a decrease in a rate or efficiency, a log-logistic relationship is assumed:

$$
y_{z,j} = \frac{1}{1 + \left( \frac{D_{z,j}}{e_{z,j}}\right)^{b}}
$$

$y_{z,j}$ is the relative response, i.e. the factor which is multiplied with a baseline value in the model.

$e_{z,j}$ is the median effective damage and $b$ is the Hill's slope. 
This formulation is not the default in most DEB-TKTD models, but generally very common in ecotoxicology.

If the stressor causes an increase in a metabolic rate (e.g. increase in maintenance costs), 
$y_{z,j}$ needs to range from 1 to $\infty$. 
This can be achieved by using the inverse of the log-logistic function, leading to the following simple relationship:

$$
y_{z,j} = 1 + \left( \frac{D_{z,j}}{e_{z,j}} \right)^{b}
$$

The parameter $e_{z,j}$ is then the damage level at which a stressor causes a doubling of a metabolic rate.

A stressor may also cause lethal effects according to the GUTS-RED-SD model (**citation needed**).

In this case, we use the cumulative hazard function to model the relationship between damage and hazard rate $h_{z}$:

$$
h_{z} = -ln \left( \frac{1}{1 + \left( \frac{D_{z,j}}{e_{z,j}}\right)^{b}} \right)
$$

Since the exponential transform $e^{-h_z dt}$ is used to convert the hazard rate to a survival probability over timespan $dt$, the dose-response of survival over damage at a given $dt$ again follows the well-known log-logistic dose-response function.

For mixture effects, we assume by default that relative responses combine multiplicatively (Independent Action, IA). 

## Agent-based model 

The ABM makes some very rudimentary assumptions about the rules which apply to an individual in a population-context. 

Reproduction is assumed to occur at fixed time intervals $\tau_R$. <br>

Aging is implemented amechanistically, where agents are assigned a maximum life span sampled from a truncated Normal distribution. <br>

The starvation-related survival probability $s_{f_X}$ decreases with decreasing $f_X$ following a linear relationship below a threshold value $f_{X,thr}$. At $f_X=1$, $s_{f_X}$ hits a non-zero value $s_{min}$, in order to avoid that individuals immediately die upon being deprived from food for any amount of time. 
The assumed relationship between $s_{f_X}$ and $f_X$ is thus given by 

$$
s_{f_X} = \begin{cases}
    (1 - smin) f_X/ f_{X,thr} + s_{min}\ \ if\ \ f_X < f_{X,thr} \\
    \ \\
    0\ \ otherwise
\end{cases}
$$

A plausible value of $s_{min}$ can be derived from the expected survival time under complete food deprivation. <br>

These rules are sufficient to obtain plausible patterns of population dynamics, but we do not claim that they are of a specific scientific relevance. <br>
In general, the rules applied starvation and reproduction can be highly species-specific and should be reviewed for each application case. <br>