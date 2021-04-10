# Chemgaloo: a numerical simulation program of complicated kinetic chemical reactions system

Description:
---
Chemgaloo is designed for simulating kinetics in chemical system.
It provides a more flexible way to define any number of reactions.

Usage:
---
1. import Chemgaloo
2. use 'chemical' class (called by chemgaloo.chemical) to define chemicals that participate in system
3. use 'reaction' class (called by chemgaloo.reaction) to define reactions that will occur in system
4. select appropriate reactor for your simulation, presently, only batch reactor is implemented, other
   kind of reactors will come soon!

Examples:
---
Example1: a two-step chain reaction 1 + 2 -> 3, 1 + 3 -> 4. k1 and k2 are known, so use example1 to
          see evolution of concentrations of four chemicals.
Example2: a two-step chain reaction 1 + 2 -> 3, 1 + 3 -> 4. k1 is known but k2 is unknown. Use detec
          -tor to check when will reaction system reach given ratio, and stop the whole reactions sy
          -stem
Example3: a two-step chain reaction 1 + 2 -> 3, 1 + 3 -> 4. k1 is known but k2 is unknown. By compar
          -ing experimental data, find the same time inteval both from experimental side and simulat
          -tion side, thus to find correct k2
