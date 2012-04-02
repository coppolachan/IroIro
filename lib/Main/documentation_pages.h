/*!

  @page configure Configuration and compilation
 

  Configuration and compilation was tested in the following systems:

  * - <a href="http://gcc.gnu.org/">GNU compiler</a> (tested with g++ version 4.6.1) (Linux and Apple), refer to section \ref simple
  * - <a href="http://www.open-mpi.org/">openMPI</a> (tested with version 1.4.4) for multicore version  (Linux), refer to section \ref openMPI
  * - <a href="http://software.intel.com/en-us/articles/c-compilers/">INTEL compiler</a> (tested with icpc version 12.1.2) (Linux), refer to section \ref INTEL_comp
  * - <a href="http://www-01.ibm.com/software/awdtools/xlcpp/">IBM XLC compiler</a> (tested with xlC version 11.1, cross platform compilation) (AIX), refer to section \ref AIX-XLC
  
  for bug reports and requests please send an email to \b cossu(AT)post(DOT)kek(DOT)jp 

  Other topics in this page:
  * - \ref improved_kernels
  * - \ref verbosity
  * - \ref GSL_lib
  * - \ref testing
  * - \ref Silent
  * - \ref autotools_problems

  @section help Options help

  The output of 

  @verbatim
  ./configure --help@endverbatim
  
  shows the various options that can be passed at configuration time.
  
  Anyway, in most cases, the configuration command line is very simple as explained in the following sections.


  @section simple Simple compilation
 
  
  For a single core run just compile with
  
  @verbatim
  ./configure
  make @endverbatim
  
  This is valid on x86_64 Linux systems as well as apple-darwin Unix environments.

  To override the compiler flags use the environment variables in the <tt>./configure</tt> command line. See <tt>./configure --help</tt> for more informations.  


  @section AIX-XLC Compilation on AIX with XLC
  
  Use the following lines
  
  @verbatim
  ./configure 
  make@endverbatim
  
  Default optimization level is \c -O3 .

  To override compiler flags just use the \c CXXFLAGS="..." command, remembering that the following flags are necessary for correct configuration:

  @verbatim
  -q64 -qlanglvl=stdc99 -qrtti=type@endverbatim
  
  @section INTEL_comp Compilation with INTEL compiler
  
  The configure will automatically look for INTEL compiler. If it is found in the path, \c icpc will be used in compilation, no further specification is needed.
  
  Use \c CXXFLAGS="..." during configure to setup your preferred flags (this will override every default flag).

  \b Note: The compiler must be at least version 12.0.5

  @section openMPI Compilation on MPI capable general architecture
 
  Use this on a general multicore machine with openMPI installed

  @verbatim
  ./configure --enable-mpi@endverbatim

  The configuration program will look for most common mpi compilers in your system.

  Run with <tt>mpirun -np #nodes #executable</tt>
  
  @section improved_kernels Improved kernels

  Use the command 

  @verbatim
  ./configure --enable-improved@endverbatim
  
  to enable compilation of improved version of kernels (default = disabled, slower but stable and easy to debug - use this version for reference purposes).
  

  @section verbosity Verbosity control 
  
  In order to set the verbosity of executable files during runtime use the variable \c code_verbosity in configuration.
  
  Example
  @verbatim
  ./configure code_verbosity=5@endverbatim
  
  sets maximum verbosity. Allowed values from 0 to 5. Default is 1.
  
  @section GSL_lib Installing and linking GSL

  The <a href="http://www.gnu.org/software/gsl/">Gnu Scientific Library (GSL)</a> is required during compilation of Domain Wall routines. Configuration will fail if this is not found in your build system.
  
  On AIX the default development environment is in 32 bit mode, but IroIro compilation forces 64 bit mode. In this case please check that GLS library
  is installed in 64 bit version otherwise clash on libraries names could occur. 

  If you are compiling the library in a custom installation with XLC on AIX please setup the environment variable 

  @verbatim
  OBJECT_MODE=64@endverbatim
  
  (use <tt>export OBJECT_MODE=64</tt> in bash shell), or use

  @verbatim
  -q64@endverbatim  

  among compilation flags (<tt>CFLAGS="-q64"</tt>).

  For non standard installation \c configure assumes that GSL can also reside in \c ~/gsl/ directory, so that you can just create a symbolic link to your preferred build. 

  
  @section testing Testing and debugging

  @verbatim
  valgrind --leak-check=yes --log-file=valgrind.log #executable@endverbatim
  

  @section Silent Silent compilation

  By default now the compilation verbosity is set to the minimum. In order to disable the silent rules use the command 

  @verbatim
  --disable-silent-rules@endverbatim  

  at configure time. See the output of <tt>./configure --help</tt> for further informations.


  @section autotools_problems Problems with autotools

  Sometime you need to refresh your configure script.
  In the case you need to regenerate the <tt>configure</tt> script by yourself the <tt>autoreconf</tt> command is most probably not going to work if you have an older version of autotools (<2.68).

  In this case the workaround is just the following (tested on Linux with autoconf 2.63 and AIX 7.1): first delete the follwing files 

  @verbatim
  aclocal.m4
  Makefile.in
  lib/Makefile.in@endverbatim

  then execute the following commands in sequence, ignoring the eventual warning outputs
  @verbatim
  aclocal -I m4
  automake
  autoconf  @endverbatim

  Your brand new <tt>configure</tt> file should be generated without problems. 

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
  @page Using XML input

  The code uses %XML files for controlling parameters and object creation.

  Refer to the following section for detailed description:

  - \subpage creation_basic
  - \subpage creation_HMC
  - \subpage actionPage

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!

  @page creation_basic Creating %XML for %HMC runs (Basics)
  
  The %XML file starts with the declaration of type
  
  @verbatim <?xml version="1.0" encoding="UTF-8" ?>@endverbatim

  Remember that every section (element) must be opened and closed in the following way
  @verbatim
  <SectionName>
  ..... useful stuff, eventually nested sections
  </SectionName>@endverbatim

  %XML parser is case sensitive.

  The very first entry should always be a \c \<Parameters\> 
  element that is  closed at the end of the file by \c \</Parameters\> .

  This element encloses all the parameters needed to execute the program.

  The following, mandatory section is \c \<%Geometry\> that contains the \c \<%Lattice\> and \c \<%Nodes\> elements.

  \c \<%Lattice\> element contains 4 integers separated by spaces that define the global lattice in 4 dimensions.

  \c \<%Node\> element will contain 4 integers separated by spaces that describe the machine geometry in the four directions. 

  Following \c \<%Geometry\> is the \c \<%Configuration\> section.

  \c \<%Configuration\> contains a type tag and will eventually declare the initial configuration file.
 
  A typical line is the following 

  @verbatim
  <Configuration Type="Binary">my_conf.bin</Configuration>@endverbatim

  Possible types are  
  - \b Unit (start from all unit links, filename is ignored if present)
  - \b TextFile (text file containing the configuration in one single column)
  - \b Binary (binary file - binary format is the same as ILDG but no XML header by now)
  - \b JLQCDLegacy (reads the configurations generated by IBM BlueGene/L at KEK)

  After these initial descriptions the \c \<%HMC\> section starts. See page \ref creation_HMC

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!

  @page creation_HMC Creating %XML for %HMC runs (%HMC section)

  The %HMC section describes the core of the calculation.

  A typical appearance is the following one
  @verbatim
   <HMC>
    <Nsweeps>1000</Nsweeps>
    <Thermalization>100</Thermalization>
    <StartingConfig>10</StartingConfig>
    <SaveInterval>5</SaveInterval>
    <SavePrefix>HMC_</SavePrefix>

    <RandomNumberGen name="Mersenne Twister">
      <seedFile>seed_file</seedFile>
    </RandomNumberGen>
    
    <Integrator name="leapfrog_multistep">
    .... lot of info described later
    </Integrator>
   </HMC>@endverbatim

   The sections \c \<%Nsweeps\> , \c \<%Thermalizations\> and  \c \<%SaveInterval\> do exactly what one expects: declare the number of sweeps for be performed, the number of thermalization steps to be done before that and the sweeps interval between configuration storage calls.

   \c \<%Thermalizations\> and  \c \<%SaveInterval\> are not mandatory. If omitted the code will just use default values (0 and 1 respectively). Set \c SaveInterval to 0 to avoid configuration storage. 
  
   \c \<%SavePrefix\> is used (optionally) to give a prefix to the configuration names (Default, if not provided, is \c "Conf_" , such that configuration names will be like \c "Conf_10").

   \c \<%StartingConfig\> gives the starting sweep number [default = 1].

   The \c \<%RandomNumberGen\> section describes and initializes the random number generator.
   Only one choice is available at the moment: <b>Mersenne Twister</b> generator.

   Random number generator can be provided with a seed file in the same way as the example or can be initialized by the \c \<%init\> section with a sequence of integers (in decimal or hexadecimal format).

   The presence of the  \c \<%seedFile\> section supersedes the \c \<%init\> section, that is ignored if present.

   @section integrator \<Integrator\>  section

   This section declare the whole structure for the Molecular Dynamics trajectory.

   Mandatory items are shown in the example below

   @verbatim
    <Integrator name="leapfrog_multistep">
      <MDsteps>10</MDsteps>
      <step_size>0.02</step_size>
      <exp_approx>8</exp_approx>
      <step>
	<multiplier>1</multiplier>
	<Action ...>
	...
	</Action>
        <step>
	  <multiplier>1</multiplier>
	  <Action ...>
          ...
          </Action>
	</step>
      </step>
    </Integrator>@endverbatim

    The integrator name is declared by the \c name tag. Currently available names are

    - \b leapfrog_multistep (leapfrog integrator with multiple time scales)

    The three mandatory sections are \c \<%MDsteps\> , \c \<%step_size\> , \c \<%exp_approx\> that respectively declare the number of steps in the Molecular Dynamics evolution, the step size and the degree of approximation in the exponential function.

    The definition if the Hamiltonian is given by the several action terms that can be combined together. 

    @subsection leapfrog_multi The leapfrog_multistep integrator

    In the case of multistep integrator every action belongs to a step level. Step levels are declared in a nested tree. Every nested level automatically declares a new time scale. The step size of the new time scale is controlled by the \c \<%multiplier\> section (accepting integer values): a multiplier value of 2 for example divide by 2 the step size, resulting in a finer integration of that level <em>with respect to the previous one</em>.

    Multiple actions can belong to the same level. 

    The actions are discussed in the \ref actionPage 

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!

    @page actionPage Creating %XML for %HMC runs (%Action section)

    The \c \<%Action\> section declares a piece of the Hamiltonian, fermionic or gauge piece. It contains two tags, the \c type and \c name . 

    The \c type tag can have two different values with obvious meaning
    - \b Gauge 
    - \b Fermion

    Possible names for the action are (in the current implementation) 
    - \b Gauge
         - \b Wilson (Wilson type action, class ActionGaugeWilson)
	 - \b Rectangle (Rectangle type action, class ActionGaugeRect)
	 - \b Iwasaki (Iwasaki type action, specialization of ActionGaugeRect)
	 - \b Symanzik (Symanzik type action, specialization of ActionGaugeRect)
	 - \b DBW2 (DBW2 type action, specialization of ActionGaugeRect)
    - \b Fermion
         - \b TwoFlavors (Two flavors action, class Action_Nf2)
	 - \b TwoFlavorsRatio (Two flavors ratio of operators, class Action_Nf2_ratio)
	 - \b TwoFlavorsDomainWall (%Action for two flavors of Domain Wall fermions, class Action_Nf2_DomainWall)
	 - \b TwoFlavorsEvenOdd (%Action for two flavors even-odd preconditioned, class Action_Nf2_EvenOdd)

    Each one of these will be explained in the following sections.

    @section GWilson Gauge - Wilson action

    The action is \f$ S_G = \frac{\beta}{N_c} \sum P^{1 \times 1}\f$
    
    This is easily defined by a structure like:

    @verbatim
    <Action type="Gauge" name="Wilson">
      <beta>6.0</beta>
    </Action>@endverbatim
    
    Just the \c \<%beta\> section is necessary and provides the \f$ \beta \f$ value.


    @section GRect Gauge - Rectangle action

    The action is \f$ S_G = \frac{\beta}{N_c}( \sum c_1 \cdot P^{1 \times 1} + \sum  c_2 \cdot P^{1 \times 2})\f$

    Defined in a way similar to the standard Wilson action plus the two additional parameters \c c_plaq and \c c_rect

    @verbatim
    <Action type="Gauge" name="Rectangle">
      <beta>2.0</beta>
      <c_plaq>2.0</c_plaq>
      <c_rect>1.5</c_rect>
    </Action>@endverbatim
 
    @section GRectSpecial Gauge - Iwasaki, Symanzik, DBW2 actions

    These are defined in a way identical to the Wilson action since the coefficients for the several terms are automatically defined in the code, being respectively:

    - Iwasaki action: \f$c_1 =  3.648\f$, \f$c_2 = -0.331\f$
    - Symanzik action: \f$c_1 =  5/3\f$, \f$c_2 = -1/12\f$
    - DBW2 action: \f$c_1 =  12.2704\f$, \f$c_2 = -1.4088\f$
    
    So, a construction will look like the following:
 
    @verbatim
    <Action type="Gauge" name="Iwasaki">
      <beta>2.25</beta>
    </Action>@endverbatim


    @section FActionNf2 Fermion - TwoFlavors

    In principle this action accepts all kind of fermions, but in the current implementation only Wilson fermions are fine. This is a quite simple action too, no parameters except the declaration of the kernel are necessary.

    The example code is 
    @verbatim
    <Action type="Fermion" name="TwoFlavors">
      <Kernel name="DiracWilson">
        <mass>1.0</mass>
      </Kernel>
      <Solver type="Solver_CG">
        <MaxIter>1000</MaxIter>
	<Precision>10e-8</Precision>
      </Solver>
    </Action>@endverbatim

    The \c \<%Kernel\> section is mandatory and accepts a \c name tag for the %Dirac operator object. Please refer to \ref DiracOps for the DiracWilson operator.

    The \c \<%Solver\> section is mandatory as well, and declares which kind of linear solver should be used to calculate \f$ D_w^{-1}\f$. Refer to the page \ref solversPage for further details on different solvers.

    @section FActionNf2ratio Fermion - TwoFlavorsRatio

    The \c TwoFlavorsRatio action is constructed in the following way

    @verbatim
    <Action type="Fermion" name="TwoFlavorsRatio">
      <Numerator name="your Dirac operator">
      ... dirac parameters
      </Numerator>
      
      <Denominator name="another Dirac operator">
      ... dirac parameters
      </Denominator>
      
      <SolverNumerator type="your solver name for numerator">
      ... solver parameters
      </SolverNumerator>
      
      <SolverDenominator type="your solver name for denominator">
      ... solver parameters
      </SolverDenominator>
    </Action>@endverbatim
    
    The meaning of the several sections is self-explanatory. Refer to \ref DiracOps and  \ref solversPage for details on Dirac Operators and linear Solvers.
	

    @section FActionDWF Fermion - TwoFlavorsDomainWall

    A simple example is better to explain how to construct this action

    @verbatim
    <Action type="Fermion" name="TwoFlavorsDomainWall">
      <Kernel5D>
       ...
      </Kernel5D>
      <Solver type="your favourite solver name">
       ...
      </Solver>  
    </Action>@endverbatim

    The structure is very simple. The \c \<%Kernel5D\> will automatically construct a 5D Domain Wall operator with parameters described in section \ref DDWF5d .

    For possible solver names and parameters to be provided please refer to the page for \ref solversPage . 

    @section FActionEO Fermion - TwoFlavorsEvenOdd

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!    


    @page DiracOps %Dirac operators

    Several %Dirac operators are defined in the current version of the code and can be referenced with the following names (case sensitive)

    - \b DiracWilson (Wilson fermions Dirac operator, class Dirac_Wilson)
    - \b DiracWilson_EvenOdd (Wilson fermions Dirac operator with even-odd preconditioning, class Dirac_Wilson_EvenOdd)
    - \b DiracClover (Clover type Fermions, class Dirac_Clover) 
    - \b DiracOptimalDomainWall4d (4 dimensional Domain Wall operator, class Dirac_optimalDomainWall_4D - uses the 5 dimensional operator as a member)
    - \b DiracOptimalDomainWall5d (5 dimensional representation of Domain Wall operator, class Dirac_optimalDomainWall)

    @section DWilson DiracWilson operator

    This operator needs the structure:

    @verbatim
    <... name="DiracWilson">
       <mass>1.0</mass>
    </...>@endverbatim

    Dots are used because the section name requesting the operator depends on the upper level in xml file (see for examples the \ref actionPage ). 

    Only the \c \<%mass\> section is necessary. The mass is related to \f$\kappa\f$ by the relation
    \f[ \kappa = \frac{1}{2(4+m)}\f]
    
    @section DClover DiracClover operator
    
    The structure that defines this operator is 

    @verbatim
    <... name="DiracClover">
       <mass>0.1</mass>
       <Csw>1.0</Csw>
    </...>@endverbatim

    Dots are used because the section name requesting the operator depends on the upper level in xml file (see for examples the \ref actionPage ). 

    Besides the \c \<%mass\> section common to the Wilson operator the \c \<%Csw\> section defines the value of \f$c_{\rm sw}\f$


    @section DWilsonEO DiracWilson_EvenOdd operator

    @section DDWF4d DiracOptimalDomainWall4d operator

    @section DDWF5d DiracOptimalDomainWall5d operator

    The 5 dimensional Domain Wall operator is described by several parameters like in the following code snippet:
 
    @verbatim
    <... name="DiracOptimalDomainWall5d">
      <Preconditioning>NoPreconditioner</Preconditioning>
      <N5d>6</N5d>
      <wilson_mass>-1.8</wilson_mass>
      <b>2.0</b>
      <c>0.0</c>
      <mass>0.10</mass>
      <approximation name="Zolotarev">
        <lambda_min>0.1</lambda_min>
	<lambda_max>1.5</lambda_max>
      </approximation>
     </...>@endverbatim

     The \c \<%Preconditioning\> section declares which kind of preconditioners are allowed for this operator
     Choices are
     - \b NoPreconditioner (no preconditioner is applied before matrix inversion)
     - \b LUPreconditioner (LU preconditioner is applied - \b must be used in conjunction with a preconditioned linear solver, otherwise wrong results will be obtained).

     Sections \c \<%N5d\> \c \<%wilson_mass> \c \<%b\> \c \<%c\> and \c \<%mass\> provide respectively the length of the fifth dimension, the mass of Wilson Kernel operator, the \f$b\f$ and \f$c\f$ parameters and the mass \f$m\f$ like in the following equation

     \f[ insert-equation \f]

     The \c \<%approximation\> section is used to generate the weights of the 5 dimensional layers. We have two choices, the hyperbolic tangent case (giving the usual domain wall formulation) and the Zolotarev case (giving a better approximation of the sign function). Therefore the \c name tag has two variations:
     - \b Tanh (no further parameters are needed)
     - \b Zolotarev (needs the interval of approximation like in the example)



*/
///////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
  
  @page solversPage Linear solvers
  
  Solving linear equations requires the definition of a Solver object. Typically they are needed inside actions or quark propagators. Name is always specified by the tag \c name.
  
  Possible entries for the names are:
  - \b %Solver_CG (Conjugate gradient linear solver)
  - \b %Solver_CG_Precondition (Conjugate gradient linear solver with preconditioning, must be used together with a preconditioned Dirac operator otherwise no preconditioning is applied)
  - \b %Solver_BiCGStab (Bi-Conjugate gradient stabilized linear solver)
  
  @section SolvCG Solver_CG
  
  It uses the Conjugate gradient linear solver to invert matrices and requires the \c \<%MaxIter\> and \c \<%Precision\> sections like in the following example:

  @verbatim
  <... type="Solver_CG">
    <MaxIter>1000</MaxIter>
    <Precision>1e-8</Precision>
  </...>@endverbatim

  Name of solver section is dependent on caller (upper level of xml file).

  \c MaxIter and \c Precision respectively provide the maximum number of iterations for the solver (program aborts if iterations go beyond this limit) and the precision goal for the residual, the stopping condition.
  
 
  @section SolvCGPrec Solver_CG_Precondition

  This is totally analogous to the Solver_CG declaration.

  @section SolvBiCGStab Solver_BiCGStab

  Refer to the Solver_CG for details, parameters to be declared are the same.


*/
///////////////////////////////////////////////////////////////////////////////////////////////////////
/*!

  @page smearingPage Smearing Routines

  Currently the XML support is limited to just Stout smearing on a single configuration. No HMC support yet (although it is coded and working).

  We have two basic choices in smearing routines
  - \b APE for APE-like smearing of links (class Smear_APE)
  - \b Stout for stout smearing of links using a specific kernel (explained later, class Smear_Stout)

  @section APEsm APE
  
  This is quite simple and described by this equation where \f$C_\mu(x) \f$ is the new link:

  \f[  C_\mu(x) = \sum_{\nu\neq \mu}\rho_{\mu\nu}\biggl( U_\nu(x) U_\mu(x\!+\!\hat{\nu}) U_\nu^\dagger(x\!+\!\hat{\mu})
  + U^\dagger_\nu(x\!-\!\hat{\nu}) U_\mu(x\!-\!\hat{\nu})  U_\nu(x\!-\!\hat{\nu}\!+\!\hat{\mu})\biggr) \f]
  
  The XML snippet for creating the object is:
  @verbatim
  <... type="APE">
    <rho>1.0</rho>
  </...>@endverbatim

  Again, the dots are used because the object can have different names, for example see the stout call.
  The \c rho parameter can be a list of \f$\rho_{\mu\nu}\f$ values (fastest running index, \f$\nu\f$) or just one value like in the
  example, where is assumed that all coefficients (\f$ \mu \neq \nu\f$) are identical. 

  @section StoutSM Stout

  Stout links are automatically in the \f$SU(N)\f$ gauge group (current implementation only supports \f$N=3\f$).

  Assuming that \f$C_\mu(x) \f$ is the smeared link using some base smearer (like APE) the new link \f$\Omega_\mu(x)\f$ is 
  defined in the following way:

  
  \f[
  Q_\mu(x) =  \frac{i}{2}\biggl( \Omega^\dagger_\mu(x) - \Omega_\mu(x) \biggr)
  - \frac{i}{2N} \rm{tr}\biggl(\Omega^\dagger_\mu(x) -\Omega_\mu(x)\biggr)\f]
  \f[
  \Omega_\mu(x) = C_\mu(x)\ U_\mu^\dagger(x), \quad\mbox{(no summation over $\mu$)} 
  \f]

  The xml code to generate stout smearing is 

  @verbatim
  <.. type="Stout">
    <Base type="APE">
      <rho>1.0</rho>
    </Base>
  </...>
  @endverbatim


  The \c Base object defines how the  \f$C_\mu(x) \f$ link must be obtained. No other parameters are necessary.


*/
