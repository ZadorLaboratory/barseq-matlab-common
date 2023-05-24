This document provides a simple instruction on the SRVF software, organized by Qunqun Yu.

SRVF software contain 3 directories:

1. SRVFDT, the Fisher-Rao curve registration software provided by Derek Tucker. The function for curve registration is time_warping.m. The function time_warping_oneplot.m is same as time-warping.m expect that all 6 default plots appear in a single figure, and a user input colormap is required. If you only need to do Fisher-Rao curve registration, this directory contains everything.

2. curvePNSLXS, the PNS method provided by Xiaosun Lu. The function for PNS is curvepnsLXS.m.

3. Example contains a function Example.m illustrating how to use the two softwares. If you want to learn the software quick, please go to Example -> html -> Example.html, then you will get instruction on running Example -> Example.m.

Before using these softwares, please first download J. S. Marron's softwares AllMatlab7Combined.zip on http://www.unc.edu/depts/stat-or/miscellaneous/marron/Matlab7Software/ and unzip it. The softwares we need are 'General' and 'Smoothing'. Then, add 'General', 'Smoothing' and 'SRVF' to your search path. One way to do it is: Click on 'Set Path' -> 'Add with Subfolders', select 'General', 'Smoothing' and 'SRVF', then all these three and their subfolders will be added to the search path.

Now, the software is ready to use. Enjoy!


