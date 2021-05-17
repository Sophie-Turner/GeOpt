# COMP3000
<h1>Molecular Geometry Optimisation Using Evolutionary Computation</h1>
<p>Sophie Turner's final year project</p>

<h2>Vision</h2>
<p>This program is for chemists and physicists who want to estimate and view the structures of
molecules without spending hours doing calculations or using a supercomputer. Geopt is a
program which uses machine learning to predict the shapes of theoretical molecules.</p>

<h2>Supervisor</h2>
<p>Dr. David Walker</p>

<h2>User Manual</h2>
<h3>Intended use of Geopt & limitations</h3>
<p>Geopt uses an Effective Medium Theory (EMT) energy calculator to search for optimal geometric configurations of molecules.
The EMT calculator is provided by The Atomic Simulation Environment (ASE) and has been adapted by Sophie Turner for this use.
The EMT calculator is unable to work with ions or double or triple bonds. This is because it does not consider valency or any bonds between atoms, 
only the forces on atoms caused by the presence of other atoms nearby. This is why it favours some geometries which are not chemically realistic.
EMT was originally only intended for use with solid crystal structures.</p>

<h3>Software required</h3>
<p>The following software is required to use Geopt.</p>
<h4>Python 3</h4>
<p>Windows: https://www.python.org/downloads/ 
  <br/>Linux command: sudo apt-get install python3</p>
<h4>ASE</h4>
<p>Pip command: pip install --upgrade --user ase 
  <br/>Linux command: sudo apt-get install python3-ase</p>
<h4>ttkwidgets</h4>
<p>Pip command: pip install ttkwidgets 
  <br/>Linux commands: sudo add-apt-repository ppa:j-4321-i/ttkwidgets 
  <br/>sudo apt-get update
  <br/>sudo apt-get install python3-ttkwidgets</p>

<h3>Creating a molecule</h3>
<p>Choose elements by either typing the molecular formula or selecting elements from the periodic table. 
  Use correct capitalisation for element symbols.</p>

<h3>Many-molecule evolution</h3>
<p>The many-molecule algorithm is an evolutionary algorithm 
                   which creates a population of versions of the molecule 
                   with different configurations, applies mutations to 
                   generations of versions and selects the best versions 
                   based on the lowest total potential energy. It is 
                   recommended for molecules with two to twelve atoms.</p>

<h3>Per-atom exhaustive test</h3>
<p>The per-atom exhaustive test is an algorithm which moves 
              each atom in the molecule around each other atom, 
              constantly testing for the configuration with the lowest 
             total potential energy. The amount of processing required 
             is exponential to the number of atoms and it 
             is recommended for molecules with two to six atoms.</p>

<h3>Periodic boundary conditions</h3>
<p>Periodic boundary conditions are used to approximate
          a larger system from repeating the unit cell.</p>

<h3>Population size</h3>
<p>The population size is the number of versions of the 
              molecule that will exist at the same time at each step in 
              the algorithm, in each process. A larger population makes 
              the algorithm take longer to complete but may offer more 
              potential configurations.</p>

<h3>Parallel processes</h3>
<p>The number of parallel processes is how many times the entire 
            algorithm will run. It is recommended that this number is equal 
            to the number of CPU cores in your computer. This is detected 
            and set by default. If the program in unable to detect the number of CPU cores 
            in your computer it will set the default number to one, which you can change 
              manually. A higher number may offer more potential 
            configurations, but you are advised to avoid setting this to a 
            higher number than the number of CPU cores for molecules with 
            more than two atoms.</p>

<h3>Plot data size limit</h3>
<p>This is the maximum number of data points in the datasets used to 
                create plots. A smaller number will speed up execution time and 
                may also create clearer graphs.</p>

<h3>Radial mutation distribution</h3>
<p>This alters the position of atoms and moves them randomly according to
            a distribution around other atoms or themselves. The distribution
            size determines how far apart atoms can move from one another. Moving             atoms too far apart can cause the energy calculator to treat each atom             as a standalone system.
            Moving atoms too close together can result in high energies and                   repulsion.
</p>

<h3>Permutation</h3>
<p>This swaps the positions of atoms in the molecule during an evolutionary algorithm.</p>

<h3>Crossover</h3>
<p>This combines positions of atoms from two parent molecules to form a new configuration during an evolutionary algorithm.</p>

<h2>References & Other Documentation</h2>

<h3>Atomic Simulation Environment</h3>
<a href="url">https://wiki.fysik.dtu.dk/ase/</a>

<h3>Effective Medium Theory testing calculator source</h3>
<a href="url">https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/emt.html#EMT</a>

<h3>Authors of the research and creators of the above packages</h3>
<p>Ask Hjorth Larsen, Jens Jørgen Mortensen, Jakob Blomqvist,
 Ivano E. Castelli, Rune Christensen, Marcin Dułak, Jesper Friis,
 Michael N. Groves, Bjørk Hammer, Cory Hargus, Eric D. Hermes,
 Paul C. Jennings, Peter Bjerre Jensen, James Kermode, John R. Kitchin,
 Esben Leonhard Kolsbjerg, Joseph Kubal, Kristen Kaasbjerg,
 Steen Lysgaard, Jón Bergmann Maronsson, Tristan Maxson, Thomas Olsen,
 Lars Pastewka, Andrew Peterson, Carsten Rostgaard, Jakob Schiøtz,
 Ole Schütt, Mikkel Strange, Kristian S. Thygesen, Tejs Vegge,
 Lasse Vilhelmsen, Michael Walter, Zhenhua Zeng, Karsten Wedel Jacobsen.
 The Atomic Simulation Environment—A Python library for working with atoms
 J. Phys.: Condens. Matter Vol. 29 273002, 2017</p>
 
 <h3>Geopt planner</h3>
 <a href="url">https://tasks.office.com/live.plymouth.ac.uk/en-GB/Home/Planner/#/plantaskboard?groupId=94cc8cbf-c90e-473e-a3f2-7d7d5dee52d9&planId=VmtMkciqc0GG6F--BwFrqpYAESf1</a>
 
 <h3>YouTube demonstration (January 2021)</h3>
 <a href="url">https://youtu.be/ylv4J85m95Y</a>
 
 <h3>YouTube demonstration (February 2021)</h3>
 <a href="url">https://youtu.be/llvHbYyEO6Q</a> 
 
 <h3>Fortnightly blog</h3>
 <a href="url">https://dle.plymouth.ac.uk/mod/oublog/view.php?id=922745</a> 
