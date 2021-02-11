# COMP3000
<h1>Molecular Geometry Optimisation Using Machine Learning</h1>
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

<h3>Creating a molecule</h3>
<p>Choose elements by either typing the molecular formula or selecting elements from the periodic table. 
  Use correct capitalisation for element symbols.</p>

<h3>Many-molecule Evolution</h3>
<p>The many-molecule algorithm is an evolutionary algorithm 
                   which creates a population of versions of the molecule 
                   with different configurations, applies mutations to 
                   generations of versions and selects the best versions 
                   based on the lowest total potential energy. It is 
                   recommended for molecules with three to twelve atoms.</p>

<h3>Per-atom Exhaustive Test</h3>
<p>The per-atom exhaustive test is an algorithm which moves 
              each atom in the molecule around each other atom, 
              constantly testing for the configuration with the lowest 
             total potential energy. The amount of processing required 
             is proportional to the number of atoms factorial and it 
              is recommended for molecules with two to four atoms.</p>

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

<h3></h3>
<p></p>

<h2>References</h2>
<p></p>
