[**Introduction**](https://ivaleram.github.io/GLFM/) | [**Functions**](doc_functions.html) | [**Data Structures**](doc_struct.html) | [**Demos**](demos.html) | [**FAQ**](FAQ_errors.html)

Frequently Asked Questions
---------------------------

**Which types of data can the GLFM handle?**

This model is specially suitable to deal with very heterogeneous data. There are currently five types of data defined:

* 'g': continuous real-valued
* 'p': continuous positive real-valued
* 'n': discrete count data
* 'c': discrete categorical
* 'o': discrete ordinal

Documented Errors in MATLAB
----------------------------

**ERROR WITH MEX: GSL NOT FOUND (OS X)**

> mex  -lgsl -lgmp -lgslcblas IBPsampler.cpp
> Building with 'Xcode Clang++'.
> Error using mex
> In file included from
> GLFM/src/Ccode/wrapper_matlab/IBPsampler.cpp:1:
> GLFM/src/Ccode/wrapper_matlab/IBPsampler.h:6:10:
> fatal error: 'gsl/gsl_sf_exp.h' file not found
> #include <gsl/gsl_sf_exp.h>
>          ^
>          1 error generated.

To solve it, simply specify library location (where to look for the include) with -I flag

    mex  -lgsl -I/usr/local/include -lgmp -lgslcblas IBPsampler.cpp



Documented Errors in PYTHON
----------------------------

**ERROR ASCII ENCODING INSTEAD OF UTF8**

In case you get the following error at compilation time in Python:

> UnicodeDecodeError: 'ascii' codec can't decode byte 0xcc in position 32:
> ordinal not in range(128)

This error happens when the character encoding in your python version is not utf8 (the typical encoding).

To fix it, you might try to set up a new encoding as follows:

    import sys
    sys.setdefaultencoding('utf8')

More information available here:
<http://stackoverflow.com/questions/21129020/how-to-fix-unicodedecodeerror-ascii-codec-cant-decode-byte>

The solution stated above only holds for the current session. To make it permanent, you might need to change *site.py*, the file where encoding is defined. You might replace the line:

> encoding = "ascii" # Default value set by _PyUnicode_Init()

by

> encoding = "utf8"

------------------------------

**ERROR CYTHON.DISTUTILS**

The following error occurs at compilation time when Cython has not been properly installed:

> Traceback (most recent call last):
>  File "setup.py", line 2, in <module>
>      from Cython.Distutils import Extension
>      ImportError: No module named Cython.Distutils

Please install again Cython, following our [Installation Instructions](README.html)

--------------------------

**ERROR PROXY**

> conda update --all
> Fetching package metadata ...
> 
> CondaHTTPError: HTTP None None for url <None>
> Elapsed: None

An HTTP error occurred when trying to retrieve this URL.
SSLError(SSLError("Can't connect to HTTPS URL because the SSL module is not
available.",),)

This error happens when proxy settings are ill-defined. This is a more general error that affects any other program as well. You might want to review your proxy settings. In OS X, you can check:

    System Preferences > Network > Advanced > Proxies

Further information here:
<https://github.com/ContinuumIO/anaconda-issues/issues/1326>.

