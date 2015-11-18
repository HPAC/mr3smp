# Overview #

The project provides code for a _multi-threaded_ variant of the algorithm of _Multiple Relatively Robust Representations_. In most cases it is the fastest solver for dense symmetric and Hermitian eigenproblems. The code is free to use and edit and therefore can be used by users as well as library developers to integrate the code or using it as start for there own implementation.

The MRRR algorithm is a solver for the symmetric tridiagonal eigenproblem, which lies at the heart of direct methods for dense symmetric and Hermitian eigenproblems. For convenience we therefore include, besides the tridiagonal solver, routines for the following dense problems:
  * _symmetric:_ A\*x = lambda\*x, with A=A<sup>T</sup>
  * _Hermitian:_ A\*x = lambda\*x, with A=A<sup>H</sup>
  * _generalized symmetric-definite:_ A\*x = lambda\*B\*x or A\*B\*x = lambda\*x or B\*A =  lambda\*x, with A=A<sup>T</sup>, B=B<sup>T</sup> and B definite
  * _generalized Hermitian-definite:_ A\*x = lambda\*B\*x or A\*B\*x = lambda\*x or B\*A = lambda\*x, with A=A<sup>H</sup>, B=B<sup>H</sup> and B definite

The matrices A and B can be stored either full or, although not recommended, in packed storage format. The routines for the dense problems require linking to LAPACK and an optimized BLAS. Please refer to the README file for more details.



# Download #

NOTE: Since Google Code’s download option is deprecated, only older versions can be found in the download section.

You can download the lastest version of the source code here: [mr3smp-version-1.2](http://code.google.com/p/mr3smp/downloads/detail?name=mr3smp-version-1.2.tar.gz&can=2&q=). <br>A version that only consists of C code can be found here: <a href='http://code.google.com/p/mr3smp/downloads/detail?name=mr3smp-1.2-PureC.tar.gz&can=2&q='>mr3smp-version-1.2-PureC</a>.<br>
<br>A version that uses extended precision internally can be found here:<br>
<a href='http://hpac.rwth-aachen.de/~petschow/mr3smp-extended-20140513.tar.gz'>mr3smp-extended</a>.<br>
<br>
<h1>Citing mr3smp</h1>

When you use this code, kindly reference the following paper paper:<br>
<pre><code>@article{Petschow2011:254,<br>
    author  = "Matthias Petschow and Paolo Bientinesi",<br>
    title   = "MR$^3$-SMP: A Symmetric Tridiagonal Eigensolver for Multi-Core Architectures",<br>
    journal = "Parallel Computing",<br>
    year    = 2011,<br>
    volume  = 37,<br>
    number  = 12,<br>
}<br>
</code></pre>


<h1>Mixed precision MRRR</h1>

Mixed precision routines can be used to improve the accuracy of the solver (as depicted below) at moderate cost. It can also be used to improve parallelism.<br>
Details can be found in the following <a href='http://mr3smp.googlecode.com/files/IMPROVED%20ACCURACY%20AND%20PARALLELISM%20FOR%20MRRR-BASED%20EIGENSOLVERS%20%28Technical%20Report%29.pdf'>Technical Report</a>. The figure below shows the orthogonality for test set "Application".<br>
<br>
<img src='http://mr3smp.googlecode.com/files/ort.png' />

<br>A version that uses double/extended precision can be found here: <a href='http://hpac.rwth-aachen.de/~petschow/mr3smp-extended-20140513.tar.gz'>mr3smp-extended</a>.<br>
