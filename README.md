# 3D-MCF
A 3-D Minimum Cost Flow Phase Unwrapping Algorithm Based on Closure Phase for InSAR

1.Required third party package

You should have Triangle (https://www.cs.cmu.edu/~quake/triangle.html) installed to perform our algorithm. This package is used to generate Triangulated irregular network (TIN) to connect discrete pixels.


2.Data prepation

Please read matlab/InputParameters.m for more details. Currently, the updated code only works for discrete dataset.


3.Optimization

Currently, we are using a built-in function 'linprog' in Matlab to solve our equation, which however is time and memory consuming. If you know any other solutions, please let me know (liufei.whu@gmail.com). In addition, if you are using Matlab version after 2012a (the version which we used to develop our code), you may need to add the 'optim' dir to your path as the 'linprog' function has updated and may casue some problems in the later Matlab version.






If you have used our code, please cite our paper: 

F. Liu and B. Pan, "A New 3-D Minimum Cost Flow Phase Unwrapping Algorithm Based on Closure Phase," in IEEE Transactions on Geoscience and Remote Sensing, vol. 58, no. 3, pp. 1857-1867, March 2020, doi: 10.1109/TGRS.2019.2949926.
