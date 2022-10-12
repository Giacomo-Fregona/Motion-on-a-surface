# Motion on a surface :8ball: :golf:
Project for the oral examination of the course "Modelli matematici" [Matematical Models] (UniTS). Task: choose the problem you wish, arising from any field of application of mathematics, and manage to find a reasonable solution with the teoretical and numerical tools seen during the course (mainly mathematical modelling and numerical methods for the resolution of differential equations). I choose a problem arising from physics:
- consider the graph of a (smooth) function IR<sup>2</sup> **&rarr;** IR. This can be seen as a particular type of surface!
- Consider a point on the surface. Suppose this point can move being forced to lie on the surface.
- Choose a vector representing the initial speed of the point,
- add gravity to the world,
- how will the position and velocity of the point evolve?

In order to solve the problem I modelled the situation applying physical principles and tecniques of differential geometry. The computation of the trajectory requires both symbolic and numerical calculus. In particular I applied RK4 method for the resolution of the issued differential equations.

#### How to play with the model
The code consists of a script that can simulate the motion of the point on a generical surface of the type described above. In order to experiment the model you can just uncomment one of the example of combinations surface - starting values - lenght of the timeline - graphical settings, execute the code and insert the stepsize (for the mesh on the timeline). Raccomanded vaues are among 0.1 and 0.005. You can also customize the options listed above.

#### Mathematics
If you wold like to see my way to obtain the solution, you can look at the .pdf file, which was the presentation of the model for the oral examination. I hope you can find it interesting!
