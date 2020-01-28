

# DAARC

DAARC: Data Analysis Application for Radio Continuum

DAACR provides an end-to-end solution for radio continuum data from single dish cross-scan observation.



## installation

Dependence:

- [numpy and scipy](https://www.scipy.org/scipylib/download.html)
- [ppgplot](https://www.github.com/junliu/ppgplot/)
- [astropy](https://www.astropy.org/)
- [ephem](https://pypi.org/project/ephem/)
- [wx](https://wxpython.org/)
- [pubsub](https://github.com/schollii/pypubsub)



Installation:

`python setup.py install`



## a brief introduction







## standard data calibration procedures

Here the standard procedures for data calibration are listed.

### Gaussian fitting

 The first step of data calibration is usually Gaussian fitting, which is sensible as the observing targets are mostly 'point-like' to the antenna beam. The equation applied for the fitting is a combination of multiple-peak Gaussian and a linear function, i.e.,

<img src="https://render.githubusercontent.com/render/math?math=F(x)=\sum_{i=0}^N A_i\cdot e^{4\cdot ln2\frac{-(x-X_0)^2}{H^2}} + K\cdot x %2B B\quad\quad\quad(\rm{Eq}.\,\, 1)">

This is essentially powerful to minimize the influence of source confusion (or side lobe) and baseline drifting. <img align="right" width="265" src="demo/dialog_fitting.png">

- data location:
  - `Raw`: the location of raw FITS data (click `Browse` to locate it)
  - `Fit`: the location of the Gaussian fitted ascii results (click `Browse` to locate it)

- parameter setup
- `Amplitude`: initial guess of amplitude (<img src="https://render.githubusercontent.com/render/math?math=\rm{A}"> in Eq. 1)
  - `Offset`: initial guess of offset (<img src="https://render.githubusercontent.com/render/math?math=\rm{X}_0"> in Eq. 1)
  - `HPBW`: initial guess for antenna beam width (<img src="https://render.githubusercontent.com/render/math?math=\rm{H}"> in Eq. 1)
  - `Tcal`: initial guess for baseline interception (<img src="https://render.githubusercontent.com/render/math?math=\rm{B}"> In Eq. 1), equivalent to the strength of injected noise diode in Kelven
  - `Num. of Peaks`: the number of Gaussians (<img src="https://render.githubusercontent.com/render/math?math=\rm{N}"> in Eq. 1)
  - `Baseline Order`: the order of baseline. By default a linear function is fitted to the baseline. Higher order of baseline is possible, however this is in priciple not expected for a sub-scan, and may introduce over fitting.
  - `Channel`: stokes (R, L, RL, LR) or their combinations (e.g., `0.55*R+0.45*L`)
  - `De-Noise`: factor (<img src="https://render.githubusercontent.com/render/math?math=n">) controlling the strength of data smoothing. If $n\neq 0$, the data is smoothed regressively with  piecewise cubic spline algorithm, data points over <img src="https://render.githubusercontent.com/render/math?math=n\times rms"> are filtered. The lower value of <img src="https://render.githubusercontent.com/render/math?math=n"> the stronger smoothing. It is recommended that <img src="https://render.githubusercontent.com/render/math?math=n\geq 3"> even if strong smoothing is needed.
  - `Data Cut`: Restricts the range of the Gaussian fit in scanning direction on the left and right. By default both values are 0.1, which means that 10% of data are reduced in both sides.
  - `window`: factor (n) defines the location of the first null of the antenna beam pattern. The location is calculated as <img src="https://render.githubusercontent.com/render/math?math=\rm{X}_0\pm \rm{n}\cdot \rm{H}">, where <img src="https://render.githubusercontent.com/render/math?math=\rm{X}_0"> and <img src="https://render.githubusercontent.com/render/math?math=\rm{H}"> are the offset and HPBW from Gassian fitting, respectively.
  - `Scan Num.` starting and ending scans to be fitted.
- buttons

  - `browse`: locates raw FITS or fit format files
  - `Cancel`: cancel the settings and close the fitting dialog

  - `Reset`: reset the fitting parameters
  - `Go`: start Gaussian fitting

  


### pointing correction

The expected amplitude is normally underestimated due to pointing erros of the antenna. The amplitude loss indueced by such pointing offests can be easily evaluated by a 2-D Gaussian function, in which the corrected amplitude is given by:

 <img src="https://render.githubusercontent.com/render/math?math=A_{||}^{corr}=A_{||}^{obs}\cdot \exp(4\cdot ln2 \cdot \frac{X_{0\perp}^2}{H_{\perp}^2})\quad\quad\quad(\rm{Eq}.\,\, 2)">

<img align="right" width="265" src="demo/dialog_pointing.png"> Where <img src="https://render.githubusercontent.com/render/math?math=\rm{A}_{||}^{corr}">: the corrected amplitude on scan direction AZ (EL)

<img src="https://render.githubusercontent.com/render/math?math=\rm{A}_{||}^{obs}">: the observed amplitude on scan direction AZ (EL), equivelent to Gaussian parameter <img src="https://render.githubusercontent.com/render/math?math=\rm{A}"> in Eq. 1

<img src="https://render.githubusercontent.com/render/math?math=\rm{X}_{0\perp}"> : the pointing offset on scan direction EL (AZ), equivelent to Gaussian parameter <img src="https://render.githubusercontent.com/render/math?math=\rm{X}_0"> in Eq. 1

<img src="https://render.githubusercontent.com/render/math?math=\rm{H}_{\perp}"> : the pointing offset on scan direction EL (AZ), equivelent to Gaussian parameter <img src="https://render.githubusercontent.com/render/math?math=\rm{H}"> in Eq. 1

- parameter setup
  - `Update Result`: whether or not to update the joined fitted files
  - `Quality Control`: whether or not to perform data quality check for a subscan
- parameters for quality control (only applies when the flag of `Quality Control` is swithed on, AKA, to 1):
  - `Beam Width`: the expected beam width of the antenna at the observing frequency. A theoretical value will be given by clicking on the `Beam_Width` botton
  - `D_Off`: pointing offset tolerance in arcsec, by deault it is set to 100, which means that a subscan with <img src="https://render.githubusercontent.com/render/math?math=|\rm{X}_0| \gt 100"> is considered as a bad scan
  - `D_Width`:
  - `Rela_Err`:
  - `D_Symmetry`:
  - `D_Avg`:
- 



### air opacity correction

![dialog_tau](demo/dialog_tau.png)



### gain-elevation correction

![dialog_gain](demo/dialog_gain.png)



### gain-time correction

![dialog_time](demo/dialog_time.png)





### absolute flux density conversion

![dialog_flux](demo/dialog_flux.png)

